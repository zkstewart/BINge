#! python3
# BINge.py
# BIN Genes for Expression analyses

# This script/program aims to provide an alternative to Corset
# for situations where a reference genome is available. However,
# the reference genome might be for a different subspecies than
# the one you're working with. Hence, it's not good enough to
# just map against the reference genome. You should map against a
# de novo transcriptome, but leverage the genomic information
# to group transcripts into genes. That's what this does.

import os, argparse, sys, pickle
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.bin_handling import generate_bin_collections, populate_bin_collections, \
    iterative_bin_self_linking, multithread_bin_splitter
from modules.fasta_handling import FastaCollection
from modules.clustering import cluster_unbinned_sequences

# Define functions
def validate_args(args):
    # Validate input file locations
    for fastaFile in args.fastaFiles:
        if not os.path.isfile(fastaFile):
            print(f'I am unable to locate the transcript FASTA file ({fastaFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for annotationFile in args.annotationFiles:
        if not os.path.isfile(annotationFile):
            print(f'I am unable to locate the genome GFF3 file ({annotationFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for gmapFile in args.gmapFiles:
        if not os.path.isfile(gmapFile):
            print(f'I am unable to locate the GMAP GFF3 file ({gmapFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate the input files are logically sound
    if len(args.annotationFiles) != 0:
        if len(args.annotationFiles) != len(args.gmapFiles):
            print("Your genome annotation and GMAP files are incompatible!")
            print(f"I'm seeing {len(args.annotationFiles)} annotation files and {len(args.gmapFiles)} GMAP files")
            print("These numbers should be the same. You need to fix this up and try again.")
            quit()
    else:
        print("Note that you are choosing to run BINge without an input annotation")
        print("That's okay if there isn't one available, but if there IS one available, " + 
              "I'd recommend that you use it.")
        print("Program will continue operation as usual.")
    
    # Validate optional BINge parameters
    if args.threads < 1:
        print("--threads should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if args.convergenceIters < 1:
        print("--convergence_iters should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if not 0.0 < args.identity <= 1.0:
        print("--identity should be given a value greater than zero, and equal to " + 
              "or less than 1")
        print("Fix this and try again.")
        quit()
    
    # Specifically handle gmapIdentity
    if args.gmapIdentity == []:
        args.gmapIdentity = [0.95 for _ in range(len(args.gmapFiles))]
    if len(args.gmapIdentity) != len(args.gmapFiles):
        print("--gmapIdentity parameter must have the same number of values as -gm")
        print("Fix this and try again.")
        quit()
    for gmapID in args.gmapIdentity:
        if not 0.0 < gmapID <= 1.0:
            print("--gmapIdentity should be given values greater than zero, and equal to " + 
                "or less than 1")
            print("Fix this and try again.")
            quit()
    
    # Validate optional MMseqs2 parameters
    if args.unbinnedClusterer in ["mmseqs-cascade", "mmseqs-linclust"]:
        if args.mmseqsDir == None:
            print(f"--mmseqs must be specified if you're using '{args.unbinnedClusterer}'")
            quit()
        if not os.path.isdir(args.mmseqsDir):
            print(f'I am unable to locate the MMseqs2 directory ({args.mmseqsDir})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        "--tmpDir is validated by the MM_DB Class"
        if args.evalue < 0:
            print("--evalue must be greater than or equal to 0")
            quit()
        if not 0 <= args.coverage <= 1.0:
            print("--coverage must be a float in the range 0.0 -> 1.0")
            quit()
        "--mode is controlled by argparse choices"
        
        # Validate "MMS-CASCADE" parameters
        if args.sensitivity in ["5.7", "7.5"]:
            args.sensitivity = float(args.sensitivity)
        else:
            args.sensitivity = int(args.sensitivity)
        
        if args.steps < 1:
            print("--steps must be greater than or equal to 1")
            quit()
    
    # Validate optional CD-HIT parameters
    if not 0.0 <= args.cdhitShortCov <= 1.0:
        print("--cdhit_shortcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not 0.0 <= args.cdhitLongCov <= 1.0:
        print("--cdhit_longcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not args.cdhitMem >= 100:
        print("--cdhit_mem should be given a value at least greater than 100 (megabytes)")
        print("Fix this and try again.")
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def find_missing_sequence_id(binCollectionList, transcriptRecords):
    '''
    Compares one or more BinCollection objects against the transcript sequences
    to see if any sequences indicated in the BinCollection do not exist in
    transcriptRecords. If any are found to be missing, the first one will be returned
    as an exemplar.
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        missingSeqID -- a string of the first sequence encountered in our BinCollections
                        that could not be found in transcriptRecords.
    '''
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            for seqID in bin.ids:
                if seqID not in transcriptRecords:
                    return seqID
    return None

def get_unbinned_sequence_ids(binCollectionList, transcriptRecords):
    '''
    Compares one or more BinCollection objects against the transcript sequences
    to see if any sequences indicated in transcriptRecords do not exist in
    any Bins.
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        unbinnedIDs -- a set containing string value for sequence IDs that were not
                       binned by BINge's main clustering process.
    '''
    binnedIDs = []
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            binnedIDs.extend(bin.ids)
    binnedIDs = set(binnedIDs)
    
    unbinnedIDs = set()
    for record in transcriptRecords:
        if record.name not in binnedIDs:
            unbinnedIDs.add(record.name)
    
    return unbinnedIDs

def get_parameters_hash(args):
    '''
    Function to receive the arguments associated with BINge.py and generate a hash
    unique to this run. This hash can be used for pickle persistance of data which
    will enable the program to resume from the point that BINge's operations finish
    and the external clustering begins.
    
    Parameters:
        args -- the Argparse object associated with the main BINge.py function.
    Returns:
        paramHash -- a sha256 hash string of the parameters
    '''
    
    "These hashes are the only ones which behaviourally influence the pre-external clustering"
    HASHING_PARAMS = ["fastaFiles", "annotationFiles", "gmapFiles", "convergenceIters", "gmapIdentity"]
    
    strForHash = ""
    for param in HASHING_PARAMS:
        strForHash += str(args.__dict__[param])
    paramHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return paramHash

## Main
def main():
    showHiddenArgs = '--help-long' in sys.argv
    
    # User input
    usageShort = """Quick notes for the use of %(prog)s: 1) For each annotation GFF3 (-ga)
    you should provide a matching GMAP GFF3 alignment (-gm). 2) One or more inputs can be
    provided with -i, which are looked at internally as a single file; all sequences in the
    -ga files must be provided with -i (i.e., sequences with the same ID as the GFF3 gene IDs
    must be given) as well as any sequences used for GMAP alignment. 3) Unbinned sequences will be
    cascade clustered using MMseqs2 by default which is recommended; you can choose Linclust
    or CD-HIT if desired. 4) Some parameters are hidden in this short help format since their
    defaults are adequate; specify --help-long to see information for those options.
    """
    
    usageLong = """%(prog)s (BIN Genes for Expression analyses) is a program which bins
    de novo-assembled transcripts together on the basis of reference genome alignments.
    This might be necessary when working with multiple subspecies that are expected to
    diverge only slightly from the reference organism. By binning like this, each subspecies
    can have its gene counts compared fairly during DGE, and some of the pitfalls of
    other approaches e.g., CD-HIT are avoided.
    ###
    Prior to running BINge, you need to run GMAP alignment of your de novo transcriptome
    against one of more reference genomes. Running GMAP with parameters '-f 2 -n 6' is 
    suggested (-f 2 is required to have GFF3 format output!).
    ###
    Note 1: For each annotation GFF3 (-ga) you should provide a matching GMAP GFF3 alignment
    file (-gm). These values should be ordered equivalently.
    ###
    Note 2: You may want to provide multiple inputs with -i, one for your transcripts
    and another for the reference sequences from each of your GFF3(s). All sequences
    indicated in your -ga and -gm files need to be locateable in the files given to -i.
    ###
    Note 3: Sequences which do not align against the genome are considered to be "unbinned".
    These will be clustered with MMseqs2 cascaded clustering by default, which is the recommended
    choice. You can use Linclust (also okay) or CD-HIT (potentially very slow) if wanted.
    ###
    Note 4: You're seeing the --help-long format of this message, which means you may want to
    configure the way clustering of unbinned sequences works. Behavioural parameters of
    the algorithms can be tuned here, but the defaults are expected to work most of the time.
    The main exception is CD-HIT's memory utilisation, which probably should be set depending
    on what you have available. tldr; change these only if you know what you're doing since
    they have been set to defaults which are likely to be optimal for most use cases.
    """
    
    p = argparse.ArgumentParser(description=usageLong if showHiddenArgs else usageShort)
    # Required
    p.add_argument("-i", dest="fastaFiles",
                   required=True,
                   nargs="+",
                   help="Input transcriptome FASTA file(s) (mRNA or CDS)")
    p.add_argument("-gm", dest="gmapFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more GMAP (GFF3) files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    # Optional - BINge
    p.add_argument("-ga", "--annot", dest="annotationFiles",
                   nargs="+",
                   required=False,
                   help="Optionally, input one or more genome annotation (GFF3) files",
                   default=[])
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many threads to run when multithreading
                   is available (default==1)""",
                   default=1)
    p.add_argument("--convergence_iters", dest="convergenceIters",
                   required=False,
                   type=int,
                   help="""Optionally, specify a maximum number of iterations allowed for
                   bin convergence to be achieved (default==5); in most cases results will
                   converge in fewer than 5 iterations, so setting a maximum acts merely as
                   a safeguard against edge cases I have no reason to believe will ever
                   happen"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=5)
    p.add_argument("--gmapIdentity", dest="gmapIdentity",
                   required=False,
                   nargs="+",
                   type=float,
                   help="""Optionally, specify the identity threshold for accepting a GMAP
                   alignment for EACH value given to -gm (default==0.95 for each file); note
                   that this value operates independently of --identity and its strictness
                   should depend on whether it's aligning against the same species genome
                   (strict), against same genus genome (less strict) or different genus
                   genome (least strict)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=[])
    # Optional - program behavioural controls
    p.add_argument("--clusterer", dest="unbinnedClusterer",
                   required=False,
                   choices=["mmseqs-cascade", "mmseqs-linclust", "cd-hit"],
                   help="""Specify which algorithm to use for clustering of unbinned sequences
                   (default=='mmseqs-cascade')""",
                   default="mmseqs-cascade")
    p.add_argument("--mmseqs", dest="mmseqsDir",
                   required=False,
                   help="""If using MMseqs2-based clustering, specify the directory containing
                   the mmseqs executable""")
    p.add_argument("--cdhit", dest="cdhitDir",
                   required=False,
                   help="""If using CD-HIT clustering, specify the directory containing
                   the cd-hit-est executable""")
    p.add_argument("--identity", dest="identity",
                   required=False,
                   type=float,
                   help="""ALL CLUSTERERS: Specify the identity threshold for clustering
                   (default==0.98); this value should be strict unless you are clustering
                   multiple species' together for a DGE analysis"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.98)
    # Optional - MMseqs2
    p.add_argument("--tmpDir", dest="tmpDir",
                   required=False,
                   help="""MMSEQS: Specify the tmpDir for MMseqs2 running; default='mms2_tmp'
                   in your current working directory"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="mms2_tmp")
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="MMSEQS: Specify the evalue threshold for clustering (default==1e-3)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=1e-3)
    p.add_argument("--coverage", dest="coverage",
                   required=False,
                   type=float,
                   help="MMSEQS: Specify the coverage ratio for clustering (default==0.4)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--mode", dest="mode",
                   required=False,
                   choices=["set-cover", "connected-component", "greedy"],
                   help="MMSEQS: Specify the clustering mode (default=='connected-component')"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="connected-component")
    p.add_argument("--sensitivity", dest="sensitivity",
                   required=False,
                   choices=["4","5","5.7","6","7","7.5"],
                   help="MMSEQS-CASCADE: Specify the sensitivity value (default==5.7)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="5.7")
    p.add_argument("--steps", dest="steps",
                   required=False,
                   type=int,
                   help="""MMSEQS-CASCADE: Specify the number of cascaded clustering steps 
                   (default==3)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=3)
    # Optional - CD-HIT
    p.add_argument("--cdhit_shortcov", dest="cdhitShortCov",
                   required=False,
                   type=float,
                   help="""CDHIT: Specify what -aS parameter to provide
                   CD-HIT (default==0.4)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--cdhit_longcov", dest="cdhitLongCov",
                   required=False,
                   type=float,
                   help="""CDHIT: Specify what -aL parameter to provide
                   CD-HIT (default==0.4)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--cdhit_mem", dest="cdhitMem",
                   required=False,
                   type=int,
                   help="""CDHIT: Specify how many megabytes of memory to
                   provide CD-HIT (default==6000)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=6000)
    # Help controller
    p.add_argument("--help-long", dest="help-long",
                   action="help",
                   help="""Show all options, including those that are not
                   recommended to be changed"""
                   if not showHiddenArgs else argparse.SUPPRESS)
    
    args = p.parse_args()
    validate_args(args)
    
    # Load indexed transcripts for quick access
    "Load this upfront since it may be a memory limitation, and causing errors early is better than late"
    transcriptRecords = FastaCollection(args.fastaFiles)
    
    # Figure out what the hash of these parameters are, and what our pickle file is called
    paramHash = get_parameters_hash(args)
    pickleFile = os.path.join(
        os.path.dirname(os.path.abspath(args.outputFileName)), 
        f".{paramHash}.binge.pkl"
    )
    
    # Either load a pickle generated by previous BINge run ...
    if os.path.isfile(pickleFile) or os.path.islink(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            binCollection = pickle.load(pickleIn)
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        print(f"{pickleFile} already exists, but is not a file?")
        print("BINge expects this to be a file which it can read, or write to.")
        print("Something weird is happening, so I will exit the program now.")
        print(f"Move whatever is at the location of '{pickleFile}' then try again.")
        quit()
    # ... or begin pre-external clustering BINge
    else:
        # Parse each GFF3 into a bin collection structure
        collectionList = generate_bin_collections(args.annotationFiles, len(args.gmapFiles)) # keep genome bins separate to multi-thread later
        
        # Parse GMAP alignments into our bin collection with multiple threads
        novelBinCollection, multiOverlaps = populate_bin_collections(collectionList, args.gmapFiles,
                                                                     args.threads, args.gmapIdentity)
        
        # Merge bins resulting from fragmented annotation models
        for i in range(len(collectionList)):
            binCollection = collectionList[i].fix_fragments(multiOverlaps[i])
            collectionList[i] = binCollection
        
        # Merge gene bins together
        binCollection = collectionList[0]
        for i in range(1, len(collectionList)):
            binCollection.merge(collectionList[i])
        
        # Check that this transcriptome file contains the reference gene models
        missingSeqID = find_missing_sequence_id([binCollection, novelBinCollection], transcriptRecords)
        if missingSeqID != None:
            print(f"ERROR: '{missingSeqID}' was binned from GMAP or the annotation GFF3, " +
                "but does not exist in your input transcriptome FASTA.")
            print("A possible error is that your transcriptome lacks the reference sequences.")
            print("This error cannot be reconciled, so the program will exit now.")
            quit()
        
        # Split bins containing overlapping (but not exon-sharing) genes e.g., nested genes
        binCollection = multithread_bin_splitter(binCollection, args.threads)
        
        # Link and merge bins across genomes / across gene copies
        """Usually linking will unify multiple genomes together, but it may detect
        bins of identical gene copies and link them together which is reasonable
        since these would confound DGE to keep separate anyway"""
        
        binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
        novelBinCollection = iterative_bin_self_linking(novelBinCollection, args.convergenceIters)
        
        # Merge bin collections together and re-link bins
        """A novel bin in one genome may just be because of an absence in the annotation.
        Such an absence may not exist in another genome, and hence it will have a gene bin.
        These should be merged together to prevent having redundant bins."""
        
        binCollection.merge(novelBinCollection)
        binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
        
        # Write pickle file for potential resuming of program
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(binCollection, pickleOut)
    
    # Write binned clusters to file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines
        numClusters = 0
        for interval in binCollection:
            bin = interval.data
            for seqID in bin.ids:
                fileOut.write(f"{numClusters+1}\t{seqID}\tbinned\n")
            numClusters += 1
    
    # Cluster remaining unbinned sequences
    unbinnedIDs = get_unbinned_sequence_ids([binCollection], transcriptRecords)
    if len(unbinnedIDs) == 1:
        unbinnedClusterDict = { 0: list(unbinnedIDs)[0] }
    elif len(unbinnedIDs) > 0:
        unbinnedClusterDict = cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args)
    else:
        unbinnedClusterDict = {} # blank to append nothing to output file
    
    # Write output of clustering to file
    with open(args.outputFileName, "a") as fileOut:
        for clusterNum, clusterIDs in unbinnedClusterDict.items():
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum+numClusters+1}\t{seqID}\tunbinned\n") # clusterType = "unbinned"
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
