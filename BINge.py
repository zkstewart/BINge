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

import os, argparse, sys, pickle, platform, subprocess, json
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.bin_handling import generate_bin_collections, populate_bin_collections, \
    multithread_bin_splitter, iterative_bin_self_linking
from modules.fasta_handling import AnnotationExtractor, FastaCollection
from modules.gmap_handling import setup_gmap_indices, auto_gmapping
from modules.clustering import cluster_unbinned_sequences
from modules.validation import validate_args, validate_fasta
from Various_scripts.Function_packages.ZS_Utility import convert_windows_to_wsl_path

HASHING_PARAMS = ["inputFiles", "genomeFiles", # These hashes are the only ones which behaviourally
                  "convergenceIters", "gmapIdentity"] # influence the pre-external clustering

# Define functions
def symlinker(src, dst):
    '''
    Helper function to symlink files depending on host OS. For Linux, nothing special
    needs to occur. For Windows, developer mode is required. This function just acts as
    an exception handler.
    '''
    if platform.system() != 'Windows':
        os.symlink(src, dst)
    else:
        try:
            os.symlink(src, dst)
        except:
            print("os.symlink is not working on your Windows computer.")
            print("This means developer mode is probably not activated.")
            print("Google 'enable windows developer mode' to see how to do this.")
            print("Until then, this program will exit since symlinks are required.")
            quit()

def check_file_exists(fileLocation):
    '''
    Helper function to check if a symbolic link/file exists. Using os.path.exists() doesn't
    show an existing symlink on Windows if it was made through WSL. So we need to
    think outside the box to check if it exists or not.
    '''
    if platform.system() != 'Windows':
        return os.path.exists(fileLocation)
    else:
        fileLocation = convert_windows_to_wsl_path(fileLocation)
        cmd = ["wsl", "~", "-e", "ls", fileLocation]
        run_file_exists = subprocess.Popen(cmd, shell = True,
                                           stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        existsout, existserr = run_file_exists.communicate()
        if "No such file or directory" in existserr.decode("utf-8"):
            return False
        elif existsout.decode("utf-8").rstrip("\r\n ") == fileLocation:
            return True
        else:
            raise Exception(("symlink_exists encountered an unhandled situation; have a look " +
                             "at the stdout '" + existsout.decode("utf-8") + "' and stderr '" + 
                             existserr.decode("utf-8") + "' to make sense of this."))

def setup_working_directory(fileNames, genomeFiles, workingDirectory):
    '''
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        fileNames -- a list of strings pointing to FASTA or GFF3 files.
        genomeFiles -- a list of strings pointing to genome FASTAs to align against.
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
    '''
    numGFF3s = 0
    numFASTAs = 0
    numGenomes = 0
    
    # Create subdirectory for files (if not already existing)
    gff3Dir = os.path.join(workingDirectory, "gff3s")
    genomesDir = os.path.join(workingDirectory, "genomes")
    os.makedirs(gff3Dir, exist_ok=True)
    os.makedirs(genomesDir, exist_ok=True)
    
    # Link to the -i fileNames values
    for file in fileNames:
        # Handle GFF3:FASTA pairs
        if "," in file:
            gff3, fasta = file.split(",")
            numGFF3s += 1
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                print(f"-i value '{file}' value after the ',' is not a FASTA file")
                print("Make sure you specify the file order as GFF3:FASTA then try again.")
                quit()
            
            # Symlink files to GFF3s subdirectory if not aleady existing
            linkedGFF3 = os.path.join(gff3Dir, f"annotation{numGFF3s}.gff3")
            linkedFASTA = os.path.join(gff3Dir, f"genome{numGFF3s}.fasta")
            
            if not check_file_exists(linkedGFF3):
                symlinker(gff3, linkedGFF3)
            if not check_file_exists(linkedFASTA):
                symlinker(fasta, linkedFASTA)
        
        # Handle plain FASTA files
        else:
            numFASTAs += 1
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(file)
            if not isFASTA:
                print(f"-i value '{file}' is not a FASTA file")
                print("Make sure you specify the right file and/or location then try again.")
                quit()
            
            # Symlink to main working directory if not already existing
            linkedTranscriptome = os.path.join(workingDirectory, f"transcriptome{numFASTAs}.nucl")
            
            if not check_file_exists(linkedTranscriptome):
                symlinker(file, linkedTranscriptome)
    
    # Link to the -g genomeFiles values
    for file in genomeFiles:
        numGenomes += 1
        # Handle GFF3:FASTA pairs
        if "," in file:
            gff3, fasta = file.split(",")
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                print(f"-g value '{fasta}' value after the ',' is not a FASTA file")
                print("Make sure you specify the file order as GFF3:FASTA then try again.")
                quit()
            
            # Symlink files to genomes subdirectory if not aleady existing
            linkedGFF3 = os.path.join(genomesDir, f"annotation{numGenomes}.gff3")
            linkedFASTA = os.path.join(genomesDir, f"genome{numGenomes}.fasta")
            
            if not check_file_exists(linkedGFF3):
                symlinker(gff3, linkedGFF3)
            if not check_file_exists(linkedFASTA):
                symlinker(fasta, linkedFASTA)
        
        # Handle plain genome files
        else:
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(file)
            if not isFASTA:
                print(f"-g value '{file}' is not a FASTA file")
                print("Make sure you specify the right file and/or location then try again.")
                quit()
            
            # Symlink to main working directory if not already existing
            linkedFASTA = os.path.join(genomesDir, f"genome{numGenomes}.fasta")
            
            if not check_file_exists(linkedFASTA):
                symlinker(file, linkedFASTA)

def setup_param_cache(args, paramHash):
    '''
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        args -- the argparse object of BINge called through the main function.
        paramHash -- a string of the hash for these parameters to store in the param cache.
    '''
    # Get relevant parameters from args object
    workingDirectory = args.outputDirectory
    
    # Parse any existing param cache file
    paramCacheFile = os.path.join(args.outputDirectory, "param_cache.json")
    if os.path.exists(paramCacheFile):
        try:
            with open(paramCacheFile, "r") as fileIn:
                paramsDict = json.load(fileIn)
        except:
            raise Exception((f"'{paramCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it so I can reinitialise one."))
    else:
        paramsDict = {}
    
    # Add this program run to the paramsDict cache (if needed)
    paramsDict[paramHash] = {
        param : args.__dict__[param]
        for param in HASHING_PARAMS
    }
    
    # Write updated param cache to file
    with open(paramCacheFile, "w") as fileOut:
        json.dump(paramsDict, fileOut)

def setup_sequences(workingDirectory, isMicrobial=False):
    '''
    Will take the files within the 'gff3s' subdirectory of workingDirectory and
    produce sequence files from them.
    
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
        isMicrobial -- a boolean indicating whether the organism is a microbe and hence GFF3
                       has gene -> CDS features, rather than gene -> mRNA -> CDS/exon.
    '''
    # Locate subdirectory containing files
    gff3Dir = os.path.join(workingDirectory, "gff3s")
    assert os.path.isdir(gff3Dir), \
        f"setup_sequences failed because '{gff3Dir}' doesn't exist somehow?"
    
    # Locate all GFF3/genome pairs
    filePairs = []
    for file in os.listdir(gff3Dir):
        if file.endswith(".gff3"):
            assert file.startswith("annotation"), \
                f"GFF3 file in '{gff3Dir}' has a different name than expected?"
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".gff3")[0]
            suffixNum = filePrefix.split("annotation")[1]
            assert suffixNum.isdigit(), \
                f"GFF3 file in '{gff3Dir}' does not have a number suffix?"
            
            # Check that the corresponding genome file exists
            genomeFile = f"genome{suffixNum}.fasta"
            assert check_file_exists(os.path.join(gff3Dir, genomeFile)), \
                f"Expected to find file '{genomeFile}' in directory '{gff3Dir}' but couldn't?"
            
            # Store the pairing
            filePairs.append([os.path.join(gff3Dir, file), os.path.join(gff3Dir, genomeFile), suffixNum])
    
    # Parse out mRNA sequences from each GFF3/genome pair
    for gff3File, fastaFile, suffixNum in filePairs:
        sequenceFileName = os.path.join(workingDirectory, f"annotations{suffixNum}.nucl")
        
        # Generate file if it doesn't exist
        if not check_file_exists(sequenceFileName):
            with open(sequenceFileName, "w") as fileOut:
                seqGenerator = AnnotationExtractor(gff3File, fastaFile, isMicrobial)
                for mrnaID, exonSeq, cdsSeq in seqGenerator.iter_sequences():
                    fileOut.write(f">{mrnaID}\n{exonSeq}\n")

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
    strForHash = ""
    for param in HASHING_PARAMS:
        strForHash += str(args.__dict__[param])
    paramHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return paramHash

def _debug_pickler(objectToPickle, outputFileName):
    if not os.path.exists(outputFileName):
        with open(outputFileName, "wb") as pickleOut:
            pickle.dump(objectToPickle, pickleOut)

def _debug_loader(pickleFileName):
    with open(pickleFileName, "rb") as pickleIn:
        results = pickle.load(pickleIn)
    return results

## Main
def main():
    showHiddenArgs = '--help-long' in sys.argv
    
    # User input
    usageShort = """Quick notes for the use of %(prog)s:
    1) if providing GFF3 as input, only protein-coding mRNA features will be handled; input
    FASTA file(s) are similarly expected to be nucleotide CDS or mRNA sequences.
    2) when providing GFF3 as input, you must also specify the genome FASTA it's associated
    with; you can do that with the format 'file.gff3,file.fasta'
    3) provide as many genome files as are relevant to your analysis, alongside their
    annotations (if available and high quality) in the same style e.g., 
    'annotation.gff3,genome.fasta'
    4) Unbinned sequences will be cascade clustered using MMseqs2 by default which is
    recommended; you can choose Linclust or CD-HIT if desired. 4) Many parameters are hidden
    in this short help format since their defaults are adequate; specify --help-long to see
    information for those options.
    """
    
    usageLong = """%(prog)s (BIN Genes for Expression analyses) is a program which clusters
    (or "bins") de novo-assembled transcripts together with reference genome models if they
    are available. This might be necessary when working with multiple subspecies that are
    expected to diverge only slightly from the reference organism. By binning like this, each
    subspecies can have its gene counts compared fairly during DGE, and some of the pitfalls
    of other approaches e.g., CD-HIT are avoided.
    ###
    This program requires you to have GMAP and either MMseqs2 or CD-HIT available on your
    system. It will automatically make use of these softwares as part of its operations.
    ###
    In short, BINge operates by aligning all input files (-i) against all genomes (-g). In
    doing so, isoforms will cluster by their common best alignment positions, and synteny
    will be implicitly inferred based on where each species' best alignment positions are across
    the input genomes. This biologically minded clustering approach is likely to outperform
    methods which solely look at sequence identity.
    ###
    Note 1: You should try to provide a genome file for each species you've used in the -i
    argument or, minimally, provide the most closely related genome for each one. This isn't
    strictly necessary but you should try to do this wherever possible.
    ###
    Note 2: If you provide a GFF3 in the -i argument, you must indicate the genome FASTA it
    is associated with by linking the files together as a single string. For example you might
    provide inputs like '-i transcriptome.fasta annotation.gff3,reference.fasta'; as indicated,
    a standlone FASTA is okay, but a GFF3 must be paired with its sequence file via the 
    ',' character and in the order indicated i.e., annotation GFF3 then reference FASTA.
    ###
    Note 3: Similarly, for genome FASTAs given in the -g argument, you should either provide the
    genome on its own, or alongside its annotation GFF3. For example here you might provide inputs
    like '-g genome1.fasta annotation2.gff3,genome2.fasta'; by providing a GFF3 alongside the
    genome it will pre-seed the bins along this genome based on the annotation. If the annotation
    is of a reasonable standard this is expected to make BINge perform better.
    ###
    Note 3: Sequences which do not align against the genome are considered to be "unbinned".
    These will be clustered with MMseqs2 cascaded clustering by default, which is the recommended
    choice. You can use Linclust (also okay, trades some accuracy for some speed) or CD-HIT
    (potentially very slow and possibly least accurate) if wanted.
    ###
    Note 4: The --gmapIdentity parameter should be set in the range of ... TBD.
    ###
    Note 5: You're seeing the --help-long format of this message, which means you may want to
    configure the way clustering of unbinned sequences works. Behavioural parameters of several
    features can be tuned here, but the defaults are expected to work most of the time.
    The main exception is CD-HIT's memory utilisation, which probably should be set depending
    on what you have available if you choose to use CD-HIT. tldr; change these only if you know
    what you're doing since they have been set to defaults which are likely to be optimal for
    most use cases.
    """
    
    p = argparse.ArgumentParser(description=usageLong if showHiddenArgs else usageShort)
    # Required
    p.add_argument("-i", dest="inputFiles",
                   required=True,
                   nargs="+",
                   help="""Input transcriptome FASTA(s) and/or annotation GFF3(s) paired to
                   their genome file with ',' separator""")
    p.add_argument("-g", dest="genomeFiles",
                   required=True,
                   nargs="+",
                   help="Input genome FASTA(s) to align against")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory for intermediate and final results")
    # Optional - BINge
    p.add_argument("--gmapDir", dest="gmapDir",
                   required=False,
                   help="""If GMAP is not discoverable in your PATH, specify the directory
                   containing the mmseqs executable""")
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
                   type=float,
                   help="""Optionally, specify the identity threshold for accepting a GMAP
                   alignment (default==0.95); note that this value operates independently
                   of --identity and its strictness should depend on the largest evolutionary
                   distance you have between a file given to -i and a genome given to -g e.g.,
                   this should be strict for same species only alignment, less strict for
                   same genus alignment, and least strict for different genus alignments"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.95)
    # Optional - program behavioural controls
    p.add_argument("--microbial", dest="isMicrobial",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you are providing a GFF3
                   file from a bacteria, archaea, or just any organism in which the GFF3
                   does not contain mRNA and exon features; in this case, I expect the GFF3
                   feature to have 'gene' and 'CDS' features."""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=False)
    p.add_argument("--debug", dest="debug",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you want to generate detailed
                   logging information along the way to help with debugging."""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=False)
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
    paramHash = get_parameters_hash(args)
    print(f"# Hash associated with this run is {paramHash}")
    
    # Setup sequence working directory for analysis
    setup_working_directory(args.inputFiles, args.genomeFiles, args.outputDirectory)
    setup_param_cache(args, paramHash)
    
    # Figure out if we've already run BINge here before and exit if so
    outputFileName = os.path.join(args.outputDirectory, f"BINge_clustering_result.{paramHash}.tsv")
    if os.path.exists(outputFileName):
        print(f"A BINge output file already exists at '{outputFileName}'")
        print("This program will not overwrite an existing file.")
        print("If you want to resume an existing run, make sure to move/rename/delete this file first.")
        print("Fix this before trying again.")
        quit()
    
    # Extract mRNAs from any input GFF3 annotations
    setup_sequences(args.outputDirectory, args.isMicrobial)
    
    # Establish GMAP indexes
    setup_gmap_indices(args.outputDirectory, args.gmapDir, args.threads)
    
    # Perform GMAP mapping
    gmapFiles = auto_gmapping(args.outputDirectory, args.gmapDir, args.threads)
    _debug_pickler(gmapFiles, os.path.join(args.outputDirectory, f"{paramHash}.gmapFiles.pkl"))
    #gmapFiles = _debug_loader(os.path.join(args.outputDirectory, f"{paramHash}.gmapFiles.pkl"))
    
    # Figure out what our pickle file is called
    pickleFile = os.path.join(args.outputDirectory, f"{paramHash}.binge.pkl")
    
    # Either load a pickle generated by previous BINge run ...
    if os.path.isfile(pickleFile) or os.path.islink(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            binCollection = pickle.load(pickleIn)
        print("Note: an existing pickle file with a matching parameters hash is present " + 
              f"i.e., '{pickleFile}'.")
        print("I will load that in now and resume program operation.")
    
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        print(f"{pickleFile} already exists, but is not a file?")
        print("BINge expects this to be a file which it can read or write to.")
        print("Something weird is happening, so I will exit the program now.")
        print(f"Move whatever is at the location of '{pickleFile}' then try again.")
        quit()
    
    # ... or begin pre-external clustering BINge
    else:
        # Set up a bin collection structure for each genome
        collectionList = generate_bin_collections(args.outputDirectory, args.threads)
        if args.debug:
            print(f"# Generated a list with {len(collectionList)} collections")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} contains {len(_cl)} bins")
        _debug_pickler(collectionList, os.path.join(args.outputDirectory, f"{paramHash}.collectionList.setup.pkl"))
        #collectionList = _debug_loader(os.path.join(args.outputDirectory, f"{paramHash}.collectionList.setup.pkl"))
        
        # Parse GMAP alignments into our bin collection with multiple threads
        collectionList, multiOverlaps = populate_bin_collections(collectionList, gmapFiles,
                                                                 args.threads, args.gmapIdentity)
        if args.debug:
            print(f"# Populated collections based on GMAP alignments")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} now contains {len(_cl)} bins")
        _debug_pickler([collectionList, multiOverlaps], os.path.join(args.outputDirectory, f"{paramHash}.collectionList.populated.pkl"))
        #_debug_loader()
        
        # Split bins to separate non-overlapping gene models
        collectionList = multithread_bin_splitter(collectionList, args.threads)
        if args.debug:
            print(f"# Split bins based on lack of overlap")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} now contains {len(_cl)} bins")
        _debug_pickler(collectionList, os.path.join(args.outputDirectory, f"{paramHash}.collectionList.fragmentfixed.pkl"))
        #_debug_loader()
        
        # Merge gene bins together
        binCollection = collectionList[0]
        for i in range(1, len(collectionList)):
            binCollection.merge(collectionList[i])
        if args.debug:
            print(f"# Merged all the separate bins together")
            print(f"# The combined collection now contains {len(binCollection)} bins")
        _debug_pickler(binCollection, os.path.join(args.outputDirectory, f"{paramHash}.binCollection.squashed.pkl"))
        #_debug_loader()
        
        # Merge bin collections across genomes / across gene copies
        """Usually linking will unify multiple genomes together, but it may detect
        bins of identical gene copies and link them together which is reasonable
        since these would confound DGE to keep separate anyway."""
        
        binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
        if args.debug:
            print(f"# Iteratively self linked bins based on ID sharing")
            print(f"# The combined collection now contains {len(binCollection)} bins")
        
        # Write pickle file for potential resuming of program
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(binCollection, pickleOut)
    
    # Write binned clusters to file
    with open(outputFileName, "w") as fileOut:
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
    transcriptRecords = FastaCollection([
        os.path.join(args.outputDirectory, f)
        for f in os.listdir(args.outputDirectory)
        if f.endswith(".nucl")
    ])
    unbinnedIDs = get_unbinned_sequence_ids([binCollection], transcriptRecords)
    if args.debug:
        print(f"# There are {len(unbinnedIDs)} unbinned sequences for external clustering")
    
    if len(unbinnedIDs) == 1:
        unbinnedClusterDict = { 0: list(unbinnedIDs)[0] }
    elif len(unbinnedIDs) > 0:
        unbinnedClusterDict = cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args)
    else:
        unbinnedClusterDict = {} # blank to append nothing to output file
    
    # Write output of clustering to file
    with open(outputFileName, "a") as fileOut:
        for clusterNum, clusterIDs in unbinnedClusterDict.items():
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum+numClusters+1}\t{seqID}\tunbinned\n") # clusterType = "unbinned"
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
