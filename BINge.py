#!/usr/bin/env python3
# BINge.py
# BIN Genes for Expression analyses

# This script/program aims to provide an alternative to Corset
# for situations where a reference genome is available. However,
# the reference genome might be for a different subspecies than
# the one you're working with. Hence, it's not good enough to
# just map against the reference genome. You should map against a
# de novo transcriptome, but leverage the genomic information
# to group transcripts into genes. That's what this does.

import os, argparse, sys, pickle, json
import concurrent.futures
from hashlib import sha256
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.bins import BinBundle
from modules.bin_handling import generate_bin_collections, populate_bin_collections
from modules.gmap_handling import auto_gmapping
from modules.clustering import cluster_unbinned_sequences
from modules.parsing import BINge_Results
from modules.validation import validate_args, validate_init_args, \
    validate_cluster_args, validate_view_args, validate_fasta, \
    check_for_duplicates, handle_symlink_change, touch_ok
from modules.fasta_handling import FastaCollection
from modules.setup import TargetGenome, AnnotatedGenome, Transcriptome, \
    inputs_to_json, json_to_inputs
from _version import __version__

HASHING_PARAMS = ["identity", "clusterVoteThreshold"]

# Define functions
def setup_working_directory(targetGenomeFiles, annotatedGenomeFiles, txomeFiles, locations):
    '''
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        targetGenomeFiles -- a list of strings pointing to genome FASTAs to align against;
                             may occur as individual FASTA files or as 'GFF3,genome' pairs.
        annotatedGenomeFiles -- a list of strings occurring as 'GFF3,genome' pairs.
        txomeFiles -- a list of strings pointing to transcriptome FASTAs as individual mRNA/CDS
                      files or as 'mRNA,CDS,protein' trios.
        
        locations -- a Locations object with attributes for directory locations.
    Returns:
        targetGenomes -- a list of TargetGenome objects
        annotatedGenomes -- a list of AnnotatedGenome objects; list may be empty
        transcriptomes -- a list of Transcriptome objects; list may be empty
    '''
    # Create subdirectories for files (if not already existing)
    os.makedirs(locations.sequencesDir, exist_ok=True)
    os.makedirs(locations.gff3Dir, exist_ok=True)
    os.makedirs(locations.txDir, exist_ok=True)
    os.makedirs(locations.genomesDir, exist_ok=True)
    os.makedirs(locations.mappingDir, exist_ok=True)
    
    # Link to the --ig GFF3,genome values
    annotatedGenomes = []
    for index, file in enumerate(annotatedGenomeFiles):
        gff3, fasta = [ os.path.abspath(f) for f in file.split(",") ]
        annotatedGenome = AnnotatedGenome(index+1, locations, fasta, gff3)
        annotatedGenomes.append(annotatedGenome)
    
    # Link to the --ix transcript FASTA files
    transcriptomes = []
    for index, file in enumerate(txomeFiles):
        files = [ os.path.abspath(f) for f in file.split(",") ]
        if len(files) == 1:
            txome = Transcriptome(index+1, locations, mrnaFasta=files[0], cdsFasta=None, protFasta=None)
        else:
            txome = Transcriptome(index+1, locations, mrnaFasta=files[0], cdsFasta=files[1], protFasta=files[2])
        transcriptomes.append(txome)
    
    # Link to the -i targetGenomeFiles values
    targetGenomes = []
    for index, file in enumerate(targetGenomeFiles):
        files = [ os.path.abspath(f) for f in file.split(",") ]
        if len(files) == 1:
            targetGenome = TargetGenome(index+1, locations, files[0], gff3=None)
        else:
            targetGenome = TargetGenome(index+1, locations, files[1], gff3=files[0])
        targetGenome.length_index()
        
        targetGenomes.append(targetGenome)
    
    return targetGenomes, annotatedGenomes, transcriptomes

def get_unbinned_sequence_ids(bingeResults, eliminatedIDs, transcriptRecords):
    '''
    Compares one or more BinBundle objects against the transcript sequences
    to see if any sequences indicated in transcriptRecords do not exist in
    any Bins.
    
    Parameters:
        bingeResults -- a BINge_Results object with .binned dictionary
        eliminatedIDs -- a set containing strings of sequence IDs which are not to be
                         considered for clustering.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        unbinnedIDs -- a set containing string value for sequence IDs that were not
                       binned by BINge's main clustering process.
    '''
    binnedIDs = []
    for seqIDs in bingeResults.binned.values():
        binnedIDs.extend(seqIDs)
    binnedIDs = set(binnedIDs)
    binnedIDs = binnedIDs.union(eliminatedIDs)
    
    unbinnedIDs = set()
    for record in transcriptRecords:
        if record.name not in binnedIDs:
            unbinnedIDs.add(record.name)
    
    return unbinnedIDs

def get_parameters_hash(args, hashLen=20):
    '''
    Function to receive the arguments associated with BINge.py and generate a hash
    unique to this run. This hash can be used for pickle persistance of data which
    will enable the program to resume from the point that BINge's operations finish
    and the external clustering begins.
    
    Parameters:
        args -- the Argparse object associated with the main BINge.py function.
        hashLen -- an integer indicating the length of the hash to generate (default=20)
    Returns:
        paramHash -- a sha256 hash string of the parameters
    '''
    strForHash = ""
    for param in HASHING_PARAMS:
        strForHash += str(args.__dict__[param])
    paramHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return paramHash[0:hashLen]

## Main
def main():
    showLongArgs = '--help-long' in sys.argv
    if showLongArgs:
        sys.argv.append('--help') # hack to make BINge.py --help-long work
    
    def hide(description):
        return description if showLongArgs else argparse.SUPPRESS
    
    # Help messages
    initShort = """Quick notes for the use of BINge 'initialise':
    1) if providing GFF3 as input (--ig), only protein-coding mRNA features will be handled; input
    FASTA file(s) given to --ix are similarly expected to be nucleotide CDS or mRNA sequences.
    2) when providing GFF3 as input to -i or --ig, you must also specify the genome FASTA it's associated
    with; you can do that with the format 'file.gff3,file.fasta' and in that order.
    3) provide as many genome files as are relevant to your analysis alongside their
    annotations (if available and high quality) in the same style e.g., 
    'genome.gff3,genome.fasta'.
    """
    
    initLong = """BINge operates by aligning all input files (--ig, --ix) against all genomes (-i). In
    doing so, orthology will be implicitly inferred for sequences based on their best alignment
    position across the input genomes with overlapping sequences clustering together.
    ###
    Note 1: You should try to provide a genome file for each species you've used in the --ig/--ix
    arguments or, minimally, provide the most closely related genome for each one. This isn't
    strictly necessary but you should try to do this if possible.
    ###
    Note 2: If you provide a GFF3 in the --ig argument, you must indicate the genome FASTA it
    is associated with by linking the files together as a single string. For example you might
    provide inputs like '--ig annotation.gff3,reference.fasta'
    ###
    Note 3: If you provide transcriptome FASTAs in the --ix argument, you should provide them
    as individual mRNA/CDS sequences or as trios of mRNA/CDS/protein files with ',' separator.
    For example you might provide inputs like '--ix mRNA.fasta,CDS.fasta,protein.fasta'
    ###
    Note 4: Similarly, for genome FASTAs given in the -i argument, you should either provide the
    genome on its own, or alongside its annotation GFF3. For example here you might provide inputs
    like '-i genome1.fasta genome2.gff3,genome2.fasta'; by providing a GFF3 alongside the
    genome it will pre-seed the bins along this genome based on the annotation. If the annotation
    is of a reasonable standard this is expected to make BINge perform better.
    ###
    """
    
    clusterShort = """Quick notes for the use of BINge 'cluster':
    1) Unbinned sequences will be cascade clustered using MMseqs2 by default which is
    recommended; you can choose Linclust or CD-HIT if desired.
    2) Many parameters are hidden in this short help format since their defaults are adequate;
    specify --help-long to see information for those options.
    """
    
    clusterLong = """BINge 'cluster' follows on from 'initialise' and performs the clustering
    analysis.
    ###
    Note 1: BINge requires you to install GMAP and either MMseqs2 or CD-HIT.
    ###
    Note 2: Sequences which do not align against the genome are considered to be "unbinned".
    These will be clustered with MMseqs2 cascaded clustering by default, which is the recommended
    choice. You can use Linclust (also okay, trades some accuracy for some speed) or CD-HIT
    (potentially slow and possibly least accurate) if wanted.
    ###
    Note 3: The --identity parameter should be set in the range of 0.92 to 0.98 depending
    on the evolutionary distance between your input files and the genomes you're aligning
    against. For same species, use 0.98. For species in the same genus, use 0.95.
    If you're aligning against a different genus, consider 0.92. This setting is to make sure
    that MMseqs2 clustering of unbinned sequences occurs optimally.
    ###
    Note 4: The --clusterVoteThreshold, based on objective evidence, should be set to 0.5 or
    0.66. Using a lower value will group more sequences together into fewer clusters, and a higher
    value will group less sequences together into more clusters. A value of 0.66 is likely to be more
    biologically correct based on objective evaluation, but the difference is very marginal.
    ###
    Note 5: You're seeing the --help-long format of this message, which means you may want to
    configure the way clustering of unbinned sequences works. Behavioural parameters of several
    features can be tuned here, but the defaults are expected to work most of the time.
    The main exception is CD-HIT's memory utilisation, which probably should be set depending
    on what you have available if you choose to use CD-HIT. tldr; change these only if you know
    what you're doing since they have been set to defaults which are likely to be optimal for
    most use cases.
    """
    
    viewDescription = """BINge 'view' allows you to see the metadata of an analysis directory.
    This includes details set during 'initialise' and, if available, the results of 'cluster'.
    """
    
    mainDescription = """%(prog)s pipelines the first steps of a BINge clustering analysis.
    Set up a working directory with 'initialise' and then run 'cluster' to perform the
    clustering analysis. Use 'view' to see the metadata of an analysis directory."""
    
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("--help-long", dest="help-long",
                   action="help",
                   help="""Show all options, including those that are not
                   recommended to be changed""")
    p.add_argument("-d", dest="workingDirectory",
                   required=True,
                   help="Specify the location where the analysis is or will be performed")
    p.add_argument("-v", "--version",
                   action="version",
                   version="BINge.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=mainDescription)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="BINge.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    iparser = subparsers.add_parser("initialise",
                                    aliases=["init"],
                                    parents=[p],
                                    add_help=False,
                                    help="Initialise a working directory for the BINge pipeline",
                                    description=initLong if showLongArgs else initShort)
    iparser.set_defaults(func=imain)
    
    cparser = subparsers.add_parser("cluster",
                                    parents=[p],
                                    add_help=False,
                                    help="Run BINge clustering analysis",
                                    description=clusterLong if showLongArgs else clusterShort)
    cparser.set_defaults(func=cmain)
    
    vparser = subparsers.add_parser("view",
                                    parents=[p],
                                    add_help=False,
                                    help="View the metadata of an analysis directory",
                                    description=viewDescription)
    vparser.set_defaults(func=vmain)
    
    # Init-subparser arguments
    ## Required
    iparser.add_argument("-i", dest="targetGenomeFiles",
                         required=True,
                         nargs="+",
                         help="""Input one or more genome FASTA files to use as a clustering reference;
                         can be provided as an individual FASTA file and/or as a 'file.gff3,file.fasta' pair""")
    iparser.add_argument("--ig", dest="inputGff3Files",
                         required=False,
                         nargs="+",
                         help="""Optionally, input one or more gene model annotations as 'file.gff3,file.fasta'
                         pairs; extracted gene model sequences will be clustered but the genome FASTA file
                         will not be used as a clustering reference""",
                         default=[])
    iparser.add_argument("--ix", dest="inputTxomeFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, input one or more gene model / transcriptome sequences as
                         individual FASTA files (i.e., mRNA or CDS) and/or as trios of
                         'mrna.fasta,cds.fasta,protein.fasta'""",
                         default=[])
    ## Optional (shown)
    iparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run when multithreading
                         is available (default=1)""",
                         default=1)
    iparser.add_argument("--gmapDir", dest="gmapDir",
                         required=False,
                         help="""Optionally, if GMAP is not discoverable in your PATH, specify the
                         directory containing the 'gmap' and 'gmap_build' executables""")
    iparser.add_argument("--microbial", dest="isMicrobial",
                         required=False,
                         action="store_true",
                         help="""Optionally, specify this flag if you are inputting GFF3
                         files from bacteria, archaea, or any other organism in which the GFF3
                         does not contain mRNA and exon features; in this case, the GFF3 features
                         will occur with a parental 'gene' and child 'CDS' features with no
                         'mRNA' intermediary.""",
                         default=False)
    iparser.add_argument("--translation", dest="translationTable",
                         required=False,
                         type=int,
                         help="""Optionally, specify the NCBI translation table number of your organism(s)
                         if they have a codon table which is not The Standard Code applicable for almost
                         all eukaryotes (default=1)""",
                         default=1)
    
    # Cluster-subparser arguments
    ## Optional (shown)
    cparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run when multithreading
                         is available (default=1)""",
                         default=1)
    cparser.add_argument("--gmapDir", dest="gmapDir",
                         required=False,
                         help="""Optionally, if GMAP is not discoverable in your PATH, specify the
                         directory containing the 'gmap' and 'gmap_build' executables""")
    cparser.add_argument("--clusterer", dest="unbinnedClusterer",
                         required=False,
                         choices=["mmseqs-cascade", "mmseqs-linclust", "cd-hit"],
                         help="""Optionally, specify which algorithm to use for clustering of unbinned
                         sequences (default='mmseqs-cascade')""",
                         default="mmseqs-cascade")
    cparser.add_argument("--mmseqsDir", dest="mmseqsDir",
                         required=False,
                         help="""Optionally, if using MMseqs2-based clustering and 'mmseqs' is not
                         discoverable in your path, specify the directory containing the
                         mmseqs executable""")
    cparser.add_argument("--cdhit", dest="cdhitDir",
                         required=False,
                         help="""Optionally, if using CD-HIT clustering and 'cd-hit-est' is not
                         discoverable in your path, specify the directory containing the
                         cd-hit-est executable""")
    ## Optional (hidden)
    ### General
    cparser.add_argument("--clusterVoteThreshold", dest="clusterVoteThreshold",
                         required=False,
                         type=float,
                         help=hide("""Optionally, specify the clustering vote threshold used when
                         clustering bins based on sequence ID co-occurrence (default=0.66).
                         A higher value gives more clusters and a lower value gives fewer
                         clusters; if you want to be more biologically correct, use 0.66, but
                         if you want fewer gene clusters, use 0.5.
                         """),
                         default=0.66)
    cparser.add_argument("--identity", dest="identity",
                         required=False,
                         type=float,
                         help=hide("""ALL CLUSTERERS: Specify the identity threshold for clustering
                         (default=0.98); refer to program usage information above for help with
                         choosing an appropriate value"""),
                         default=0.98)
    cparser.add_argument("--debug", dest="debug",
                         required=False,
                         action="store_true",
                         help=hide("""Optionally, specify this flag if you want to generate detailed
                         logging information along the way to help with debugging."""),
                         default=False)
    ### MMseqs2
    cparser.add_argument("--mmseqsEvalue", dest="mmseqsEvalue",
                         required=False,
                         type=float,
                         help=hide("MMSEQS: Specify the evalue threshold for clustering (default=1e-3)"),
                         default=1e-3)
    cparser.add_argument("--mmseqsCov", dest="mmseqsCoverage",
                         required=False,
                         type=float,
                         help=hide("MMSEQS: Specify the coverage ratio for clustering (default=0.4)"),
                         default=0.4)
    cparser.add_argument("--mmseqsMode", dest="mmseqsMode",
                         required=False,
                         choices=["set-cover", "connected-component", "greedy"],
                         help=hide("MMSEQS: Specify the clustering mode (default='connected-component')"),
                         default="connected-component")
    cparser.add_argument("--mmseqsSens", dest="mmseqsSensitivity",
                         required=False,
                         choices=["4","5","5.7","6","7","7.5"],
                         help=hide("MMSEQS-CASCADE: Specify the sensitivity value (default=5.7)"),
                         default="5.7")
    cparser.add_argument("--mmseqsSteps", dest="mmseqsSteps",
                         required=False,
                         type=int,
                         help=hide("""MMSEQS-CASCADE: Specify the number of cascaded clustering steps 
                         (default=3)"""),
                         default=3)
    ### CD-HIT
    cparser.add_argument("--cdhitShortCov", dest="cdhitShortCov",
                         required=False,
                         type=float,
                         help=hide("""CDHIT: Specify what -aS parameter to provide
                         CD-HIT (default=0.4)"""),
                         default=0.4)
    cparser.add_argument("--cdhitLongCov", dest="cdhitLongCov",
                         required=False,
                         type=float,
                         help=hide("""CDHIT: Specify what -aL parameter to provide
                         CD-HIT (default=0.4)"""),
                         default=0.4)
    cparser.add_argument("--cdhitMem", dest="cdhitMem",
                         required=False,
                         type=int,
                         help=hide("""CDHIT: Specify how many megabytes of memory to
                         provide CD-HIT (default=6000)"""),
                         default=6000)
    
    # View-subparser arguments
    vparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis folder to view by its hash; if not provided,
                         the most recent analysis folder will be viewed""",
                         default="most_recent")
    
    args = subParentParser.parse_args()
    validate_args(args) # updates args.workingDirectory
    
    # Split into mode-specific functions
    if args.mode in ["initialise", "init"]:
        print("## BINge.py - Initialisation ##")
        locations = validate_init_args(args)
        imain(args, locations)
    elif args.mode == "cluster":
        print("## BINge.py - Clustering ##")
        locations = validate_cluster_args(args)
        cmain(args, locations)
    elif args.mode == "view":
        print("## BINge.py - Viewing ##")
        locations = validate_view_args(args) # sets args.runDirName
        vmain(args, locations)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def imain(args, locations):
    # Setup sequence working directory for analysis
    targetGenomes, annotatedGenomes, transcriptomes = setup_working_directory(
        args.targetGenomeFiles, args.inputGff3Files,
        args.inputTxomeFiles, locations)
    
    # Process all input types prior to GMAP alignment
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        # Handle -i files
        for targetGenome in targetGenomes:
            futures.append(
                executor.submit(targetGenome.extract_sequences, args.isMicrobial, args.translationTable)
            )
            futures.append(
                executor.submit(targetGenome.gmap_index, args.gmapDir)
            )
        # Handle --ig files
        for annotatedGenome in annotatedGenomes:
            futures.append(
                executor.submit(annotatedGenome.extract_sequences, args.isMicrobial, args.translationTable)
            )
        # Handle --ix files
        for transcriptome in transcriptomes:
            futures.append(
                executor.submit(transcriptome.extract_sequences, args.translationTable)
            )
    
    # Gather the modified argument objects
    '''
    ProcessPoolExecutor creates copies of the object in the spawned process, which means that changes
    to the .cds property within the asynchronous code do not get relayed back to the original object.
    Hence, the .extract_sequences() method returns 'self' with modifications applied. We can retrieve
    the modified object through the .result() of the process, and refer to .thisType to organise it
    into its respective list.
    '''
    targetGenomes, annotatedGenomes, transcriptomes = [], [], []
    for f in futures:
        futureResult = f.result() # raises Exceptions if any occurred
        if futureResult == None: # this is the GMAP index future
            pass
        elif futureResult.thisType == "target":
            targetGenomes.append(futureResult)
        elif futureResult.thisType == "annotated":
            annotatedGenomes.append(futureResult)
        else:
            transcriptomes.append(futureResult)
    
    # Validate that sequence duplication does not exist
    check_for_duplicates(locations.get_sequenceFiles(
        targetGenomes, annotatedGenomes,transcriptomes, "aa") # AA files are smaller than CDS or mRNA with identical IDs
    )
    
    # Perform GMAP mapping
    auto_gmapping(targetGenomes, annotatedGenomes, transcriptomes,
                  locations.mappingDir, args.gmapDir, args.threads)
    
    # Store the argument objects for use during clustering
    inputs_to_json(locations, targetGenomes, annotatedGenomes, transcriptomes, args.isMicrobial)
    
    print("Initialisation complete!")

def cmain(args, locations):
    paramHash = get_parameters_hash(args)
    print(f"# Hash associated with this 'cluster' analysis is {paramHash}")
    locations.runName = paramHash
    
    # Set up analysis directory for this run
    os.makedirs(locations.analysisDir, exist_ok=True)
    
    runDir = os.path.join(locations.analysisDir, locations.runName)
    os.makedirs(runDir, exist_ok=True)
    
    mostRecentDir = os.path.join(locations.analysisDir, "most_recent")
    if os.path.exists(mostRecentDir) or os.path.islink(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(runDir, mostRecentDir)
    
    # Load the argument objects to identify our input files
    targetGenomes, annotatedGenomes, transcriptomes, isMicrobial = json_to_inputs(locations)
    args.sequenceFiles = locations.get_sequenceFiles(targetGenomes, annotatedGenomes, transcriptomes, "cds")
    
    # Store the parameters used in this run
    with open(os.path.join(runDir, locations.parametersFile), "w") as paramOut:
        json.dump({ param: args.__dict__[param] for param in HASHING_PARAMS }, paramOut)
    
    # Figure out if we've already run BINge here before and exit if so
    outputFileName = os.path.join(runDir, locations.clusterFile)
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"A BINge output file already exists within '{runDir}'; " +
                              "this program will not overwrite an existing file. " +
                              "To resume or overwrite, move/rename/delete this file then try again.")
    
    # Figure out what our pickle file is called
    pickleFile = os.path.join(runDir, locations.pickleFile)
    
    # Either load a pickle generated by previous BINge run ...
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            clusterDict, eliminations = pickle.load(pickleIn)
        if args.debug:
            print("# An existing pickle file will be loaded to resume program operation.")
    
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        raise FileExistsError(f"'{pickleFile}' already exists within '{runDir}', but is not a file; " +
                              "move/delete/fix this to allow BINge to continue.")
    
    # ... or begin pre-external clustering BINge
    else:
        # Set up a bin collection structure for each genome
        collectionList = generate_bin_collections(targetGenomes, args.threads, isMicrobial)
        if args.debug:
            print(f"# Generated a list with {len(collectionList)} collections")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} contains {len(_cl)} bins")
        
        # Locate GMAP alignments
        gmapFiles = [
            os.path.join(locations.mappingDir, f)
            for f in os.listdir(locations.mappingDir)
            if f.endswith(".gff3")
        ]
        
        # Ensure that each GMAP alignment has a corresponding .ok file
        for gmapFile in gmapFiles:
            if not os.path.exists(gmapFile + ".ok"):
                raise FileNotFoundError(f"The GMAP file '{gmapFile}' lacks a .ok flag; rerun initialise " + 
                                        "before attempting to cluster!")
        
        # Parse GMAP alignments into our bin collection with multiple threads
        collectionList = populate_bin_collections(locations.genomesDir,
                                                  collectionList, gmapFiles,
                                                  args.threads)
        if args.debug:
            print(f"# Populated collections based on GMAP alignments")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} now contains {len(_cl)} bins")
        
        # Convert collections into a bundle
        binBundle = BinBundle.create_from_multiple_collections(collectionList)
        
        # Cluster bundles across and within genomes
        clusterDict, eliminations = binBundle.cluster_by_cooccurrence(args.clusterVoteThreshold)
        if args.debug:
            print(f"# Clustered bundles (across and within genomes) based on ID occurrence")
            print(f"# Cluster dictionary contains {len(clusterDict)} clusters")
        
        # Write pickle file for potential resuming of program
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump([clusterDict, eliminations], pickleOut)
    
    # Initialise a BINge_Results object for containing cluster results
    bingeResults = BINge_Results()
    bingeResults.binned = clusterDict
    
    # Cluster remaining unbinned sequences
    transcriptRecords = FastaCollection(args.sequenceFiles)
    unbinnedIDs = get_unbinned_sequence_ids(bingeResults, eliminations, transcriptRecords)
    if args.debug:
        print(f"# There are {len(unbinnedIDs)} unbinned sequences for external clustering")
    
    if len(unbinnedIDs) == 1:
        unbinnedClusterDict = { 0: list(unbinnedIDs) }
    elif len(unbinnedIDs) > 0:
        unbinnedClusterDict = cluster_unbinned_sequences(unbinnedIDs, transcriptRecords,
                                                         args, locations.tmpDir, runDir)
    else:
        unbinnedClusterDict = {} # blank to signal an absence of results
    
    # Store unbinned clusters in the BINge_results object
    bingeResults.unbinned = bingeResults.update_unbinned_ids(unbinnedClusterDict)
    
    # Write output of clustering to file
    bingeResults.write(outputFileName, clusterTypes="all")
    touch_ok(outputFileName)
    
    print("Clustering complete!")

def vmain(args, locations):
    # Load the parameters used in the analysis
    with open(os.path.join(args.runDir, locations.parametersFile), "r") as paramsIn:
        params = json.load(paramsIn)
    
    print(f"# Parameters for '{args.runDir}':")
    for key, value in params.items():
        print(f"{key}: {value}")
    print()
    
    # Derive the files used in the analysis
    targetGenomes, annotatedGenomes, transcriptomes, isMicrobial = json_to_inputs(locations)
    
    print("## File inputs:")
    print("# -i reference targets:")
    genomeLinks = [ x.fasta for x in targetGenomes ]
    annotLinks = [ x.gff3 for x in targetGenomes ]
    
    genomeOrigins = [ str(Path(x).resolve()) for x in genomeLinks ]
    annotOrigins = [ str(Path(x).resolve()) if x is not None else None for x in annotLinks ]
    for gLink, gOrigin, aLink, aOrigin in zip(genomeLinks, genomeOrigins, annotLinks, annotOrigins):
        print(f"{gLink} -> {gOrigin}")
        if aOrigin is not None:
            print(f"    pre-seeded by: {aLink} -> {aOrigin}")
    
    print()
    
    print("# --ig annotated gene models:")
    gff3Links = [ x.gff3 for x in annotatedGenomes ]
    fastaLinks = [ x.fasta for x in annotatedGenomes ]
    
    gff3Origins = [ str(Path(x).resolve()) for x in gff3Links ]
    fastaOrigins = [ str(Path(x).resolve()) for x in fastaLinks ]
    if len(gff3Links) == 0:
        print("None")
    else:
        for gLink, gOrigin, fLink, fOrigin in zip(gff3Links, gff3Origins, fastaLinks, fastaOrigins):
            print(f"{gLink} -> {gOrigin}")
            print(f"    associated FASTA: {fLink} -> {fOrigin}")
    print()
    
    print("# --ix transcript FASTAs:")
    mrnaLinks = [ x.mrna for x in transcriptomes ]
    cdsLinks = [ x.cds for x in transcriptomes ]
    aaLinks = [ x.aa for x in transcriptomes ]
    
    mrnaOrigins = [ str(Path(x).resolve()) for x in mrnaLinks ]
    cdsOrigins = [ str(Path(x).resolve()) for x in cdsLinks ]
    aaOrigins = [ str(Path(x).resolve()) for x in aaLinks ]
    
    if len(mrnaLinks) == 0:
        print("None")
    else:
        ongoingCount = 0
        for mrnaLink, mrnaOrigin, cdsLink, cdsOrigin, aaLink, aaOrigin in \
        zip(mrnaLinks, mrnaOrigins, cdsLinks, cdsOrigins, aaLinks, aaOrigins):
            print(f"mRNA: {mrnaLink} -> {mrnaOrigin}")
            
            if cdsLink == cdsOrigin:
                print(f"CDS (predicted): {cdsLink}")
            else:
                print(f"CDS: {cdsLink} -> {cdsOrigin}")
            
            if aaLink == aaOrigin:
                print(f"AA (predicted): {aaLink}")
            else:
                print(f"AA: {aaLink} -> {aaOrigin}")
            
            ongoingCount += 1
            if ongoingCount < len(mrnaLinks):
                print()
    print()
    
    # Parse the clustering results (if available)
    print("# Clustering result:")
    clusterFile = os.path.join(args.runDir, locations.clusterFile)
    if os.path.exists(clusterFile):
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        seqsBinned = sum([len(seqIDs) for seqIDs in bingeResults.binned.values()])
        seqsUnbinned = sum([len(seqIDs) for seqIDs in bingeResults.unbinned.values()])
        
        print(f"# Number of binned clusters: {len(bingeResults.binned)}")
        print(f"# Number of sequences in binned clusters: {seqsBinned}")
        print(f"# Number of unbinned clusters: {len(bingeResults.unbinned)}")
        print(f"# Number of sequences in unbinned clusters: {seqsUnbinned}")
    else:
        print("No clustering result available")
    
    print("Viewing complete!")

if __name__ == "__main__":
    main()
