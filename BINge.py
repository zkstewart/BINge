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

import os, argparse, sys, pickle, json
from hashlib import sha256
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.bins import BinBundle
from modules.bin_handling import generate_bin_collections, populate_bin_collections
from modules.gmap_handling import setup_gmap_indices, auto_gmapping
from modules.gff3_handling import extract_annotations_from_gff3
from modules.clustering import cluster_unbinned_sequences
from modules.parsing import BINge_Results
from modules.validation import initialise_working_directory, validate_init_args, \
    validate_cluster_args, validate_view_args, validate_fasta, \
    check_for_duplicates, handle_symlink_change, touch_ok
from modules.fasta_handling import FastaCollection, \
    generate_sequence_length_index, process_transcripts
from _version import __version__

HASHING_PARAMS = ["identity", "clusterVoteThreshold"]

# Define functions
def setup_working_directory(gff3Files, txomeFiles, targetGenomeFiles, locations):
    '''
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        gff3Files -- a list of strings pointing to GFF3,genome files.
        txomeFIles -- a list of strings pointing to transcriptome FASTAs as individual mRNA
                      files or as mRNA,CDS,protein files.
        targetGenomeFiles -- a list of strings pointing to genome FASTAs to align against.
        locations -- a Locations object with attributes for directory locations.
    '''
    # Create subdirectory for files (if not already existing)
    os.makedirs(locations.sequencesDir, exist_ok=True)
    os.makedirs(locations.gff3Dir, exist_ok=True)
    os.makedirs(locations.txDir, exist_ok=True)
    os.makedirs(locations.genomesDir, exist_ok=True)
    
    # Link to the --ig GFF3,genome values
    numIG = 0    
    for file in gff3Files:
        numIG += 1
        
        # Extract values from pair
        gff3, fasta = [os.path.abspath(f) for f in file.split(",")]
        
        # Check that FASTA is a FASTA
        isFASTA = validate_fasta(fasta)
        if not isFASTA:
            raise ValueError(f"--ig value '{file}' after the ',' is not a FASTA file")
        
        # Symlink files to GFF3s subdirectory if not aleady existing
        linkedGFF3 = os.path.join(locations.gff3Dir, f"annotations{numIG}.gff3")
        if os.path.exists(linkedGFF3):
            handle_symlink_change(linkedGFF3, gff3)
        else:
            os.symlink(gff3, linkedGFF3)
        
        linkedFASTA = os.path.join(locations.gff3Dir, f"annotations{numIG}.fasta")
        if os.path.exists(linkedFASTA):
            handle_symlink_change(linkedFASTA, fasta)
        else:
            os.symlink(fasta, linkedFASTA)
    
    # Link to the --ix transcript FASTA files
    numTX = 0
    for file in txomeFiles:
        numTX += 1
        
        # Extract values if triplicate
        files = [os.path.abspath(f) for f in file.split(",")]
        
        # Check that FASTA is a FASTA
        for i, f in enumerate(files):
            isFASTA = validate_fasta(f)
            if not isFASTA:
                raise ValueError(f"--ix value '{f}' is not a FASTA file")
            
            # Symlink to main working directory if not already existing
            suffix = "mrna" if i == 0 else "cds" if i == 1 else "aa"
            linkedTranscriptome = os.path.join(locations.txDir, f"transcriptome{numTX}.{suffix}")
            if os.path.exists(linkedTranscriptome):
                handle_symlink_change(linkedTranscriptome, f)
            else:
                os.symlink(f, linkedTranscriptome)
    
    # Link to the -i targetGenomeFiles values
    numGenomes = 0
    for file in targetGenomeFiles:
        numGenomes += 1
        # Handle GFF3:FASTA pairs
        if "," in file:
            gff3, fasta = [os.path.abspath(f) for f in file.split(",")]
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                raise ValueError(f"-i value '{fasta}' value after the ',' is not a FASTA file; " + 
                                 "Make sure you specify the file order as 'GFF3,FASTA' then try again.")
            
            # Symlink files to genomes subdirectory if not aleady existing
            linkedGFF3 = os.path.join(locations.genomesDir, f"genome{numGenomes}.gff3")
            if os.path.exists(linkedGFF3):
                handle_symlink_change(linkedGFF3, gff3)
            else:
                os.symlink(gff3, linkedGFF3)
            
            linkedFASTA = os.path.join(locations.genomesDir, f"genome{numGenomes}.fasta")
            if os.path.exists(linkedFASTA):
                handle_symlink_change(linkedFASTA, fasta)
            else:
                os.symlink(fasta, linkedFASTA)
        
        # Handle plain genome files
        else:
            fasta = os.path.abspath(file)
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                raise ValueError(f"-i value '{fasta}' is not a FASTA file")
            
            # Symlink to main working directory if not already existing
            linkedFASTA = os.path.join(locations.genomesDir, f"genome{numGenomes}.fasta")
            
            if os.path.exists(linkedFASTA):
                handle_symlink_change(linkedFASTA, fasta)
            else:
                os.symlink(fasta, linkedFASTA)
        
        # Index the genome's contig lengths if not already done
        if not os.path.exists(linkedFASTA) or not os.path.exists(f"{linkedFASTA}.ok"):
            generate_sequence_length_index(linkedFASTA)
            touch_ok(linkedFASTA)

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
        hashLen -- an integer indicating the length of the hash to generate (default==20)
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
    as individual mRNA/CDS sequences or as triplicates of mRNA/CDS/protein files with ',' separator.
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
                         help="""Input genome FASTA(s) to align against as individual files or
                         as 'file.gff3,file.fasta' pairs""")
    iparser.add_argument("--ig", dest="inputGff3Files",
                         required=False,
                         nargs="+",
                         help="""Optionally, input annotation GFF3(s) paired to their genome file
                         with ',' separator""",
                         default=[])
    iparser.add_argument("--ix", dest="inputTxomeFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, input transcriptome FASTA(s) as individual files or as triplicates
                         of mRNA/CDS/protein files with ',' separator""",
                         default=[])
    ## Optional (shown)
    iparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run when multithreading
                         is available (default==1)""",
                         default=1)
    iparser.add_argument("--gmapDir", dest="gmapDir",
                         required=False,
                         help="""If GMAP is not discoverable in your PATH, specify the directory
                         containing the 'gmap' and 'gmap_build' executables""")
    iparser.add_argument("--microbial", dest="isMicrobial",
                         required=False,
                         action="store_true",
                         help="""Optionally provide this argument if you are providing a GFF3
                         file from a bacteria, archaea, or just any organism in which the GFF3
                         does not contain mRNA and exon features; in this case, I expect the GFF3
                         feature to have 'gene' and 'CDS' features.""",
                         default=False)
    iparser.add_argument("--translation", dest="translationTable",
                         required=False,
                         type=int,
                         help="""Optionally provide this argument if you have provided individual files in
                         --ix which will be translated; default is 1 (standard code); indicate the NCBI
                         translation table number if your organisms have a different genetic code""",
                         default=1)
    
    # Cluster-subparser arguments
    ## Optional (shown)
    cparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run when multithreading
                         is available (default==1)""",
                         default=1)
    cparser.add_argument("--gmapDir", dest="gmapDir",
                         required=False,
                         help="""If GMAP is not discoverable in your PATH, specify the directory
                         containing the 'gmap' and 'gmap_build' executables""")
    cparser.add_argument("--clusterer", dest="unbinnedClusterer",
                         required=False,
                         choices=["mmseqs-cascade", "mmseqs-linclust", "cd-hit"],
                         help="""Specify which algorithm to use for clustering of unbinned sequences
                         (default=='mmseqs-cascade')""",
                         default="mmseqs-cascade")
    cparser.add_argument("--mmseqsDir", dest="mmseqsDir",
                         required=False,
                         help="""If using MMseqs2-based clustering and 'mmseqs' is not discoverable
                         in your path, specify the directory containing the mmseqs executable""")
    cparser.add_argument("--cdhit", dest="cdhitDir",
                         required=False,
                         help="""If using CD-HIT clustering and 'cd-hit-est' is not discoverable
                         in your path, specify the directory containing the cd-hit-est executable""")
    ## Optional (hidden)
    ### General
    cparser.add_argument("--clusterVoteThreshold", dest="clusterVoteThreshold",
                         required=False,
                         type=float,
                         help=hide("""Optionally, specify the clustering vote threshold used when
                         clustering bins based on sequence ID co-occurrence (default == 0.66).
                         A higher value gives more clusters and a lower value gives fewer
                         clusters; if you want to be more biologically correct, use 0.66, but
                         if you want fewer gene clusters, use 0.5.
                         """),
                         default=0.66)
    cparser.add_argument("--microbial", dest="isMicrobial",
                         required=False,
                         action="store_true",
                         help="""Optionally provide this argument if you are providing a GFF3
                         file from a bacteria, archaea, or just any organism in which the GFF3
                         does not contain mRNA and exon features; in this case, I expect the GFF3
                         feature to have 'gene' and 'CDS' features.""",
                         default=False)
    cparser.add_argument("--identity", dest="identity",
                         required=False,
                         type=float,
                         help=hide("""ALL CLUSTERERS: Specify the identity threshold for clustering
                         (default==0.98); refer to --help-long information for help with choosing an
                         appropriate value"""),
                         default=0.98)
    cparser.add_argument("--debug", dest="debug",
                         required=False,
                         action="store_true",
                         help=hide("""Optionally provide this argument if you want to generate detailed
                         logging information along the way to help with debugging."""),
                         default=False)
    ### MMseqs2
    cparser.add_argument("--mmseqsEvalue", dest="mmseqsEvalue",
                         required=False,
                         type=float,
                         help=hide("MMSEQS: Specify the evalue threshold for clustering (default==1e-3)"),
                         default=1e-3)
    cparser.add_argument("--mmseqsCov", dest="mmseqsCoverage",
                         required=False,
                         type=float,
                         help=hide("MMSEQS: Specify the coverage ratio for clustering (default==0.4)"),
                         default=0.4)
    cparser.add_argument("--mmseqsMode", dest="mmseqsMode",
                         required=False,
                         choices=["set-cover", "connected-component", "greedy"],
                         help=hide("MMSEQS: Specify the clustering mode (default=='connected-component')"),
                         default="connected-component")
    cparser.add_argument("--mmseqsSens", dest="mmseqsSensitivity",
                         required=False,
                         choices=["4","5","5.7","6","7","7.5"],
                         help=hide("MMSEQS-CASCADE: Specify the sensitivity value (default==5.7)"),
                         default="5.7")
    cparser.add_argument("--mmseqsSteps", dest="mmseqsSteps",
                         required=False,
                         type=int,
                         help=hide("""MMSEQS-CASCADE: Specify the number of cascaded clustering steps 
                         (default==3)"""),
                         default=3)
    ### CD-HIT
    cparser.add_argument("--cdhitShortCov", dest="cdhitShortCov",
                         required=False,
                         type=float,
                         help=hide("""CDHIT: Specify what -aS parameter to provide
                         CD-HIT (default==0.4)"""),
                         default=0.4)
    cparser.add_argument("--cdhitLongCov", dest="cdhitLongCov",
                         required=False,
                         type=float,
                         help=hide("""CDHIT: Specify what -aL parameter to provide
                         CD-HIT (default==0.4)"""),
                         default=0.4)
    cparser.add_argument("--cdhitMem", dest="cdhitMem",
                         required=False,
                         type=int,
                         help=hide("""CDHIT: Specify how many megabytes of memory to
                         provide CD-HIT (default==6000)"""),
                         default=6000)
    
    # View-subparser arguments
    vparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis folder to view by its hash; if not provided,
                         the most recent analysis folder will be viewed""",
                         default="most_recent")
    
    args = subParentParser.parse_args()
    initialise_working_directory(args.workingDirectory)
    
    # Split into mode-specific functions
    if args.mode in ["initialise", "init"]:
        print("## BINge.py - Initialisation ##")
        locations = validate_init_args(args)
        imain(args, locations)
    elif args.mode == "cluster":
        print("## BINge.py - Clustering ##")
        locations = validate_cluster_args(args) # sets args.sequenceFiles
        cmain(args, locations)
    elif args.mode == "view":
        print("## BINge.py - Viewing ##")
        locations = validate_view_args(args) # sets args.runDirName
        vmain(args, locations)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def imain(args, locations):
    # Setup sequence working directory for analysis
    setup_working_directory(args.inputGff3Files, args.inputTxomeFiles,
                            args.targetGenomeFiles, locations)
    
    # Extract sequences from any -i files
    extract_annotations_from_gff3(locations.gff3Dir, locations.sequencesDir,
                                  "annotations", args.threads,
                                  args.isMicrobial, args.translationTable)
    
    # Extract sequences from any --ig files
    extract_annotations_from_gff3(locations.genomesDir, locations.genomesDir,
                                  "genome", args.threads,
                                  args.isMicrobial, args.translationTable,)
    
    # Extract CDS/proteins from any --ix files
    process_transcripts(locations, args.threads, args.translationTable)
    
    # Validate that sequence duplication does not exist
    check_for_duplicates(locations.sequencesDir, ".aa") # AA files are smaller than CDS or mRNA with identical IDs
    
    # Establish GMAP indexes
    setup_gmap_indices(locations, args.gmapDir, args.threads)
    
    # Perform GMAP mapping
    auto_gmapping(locations, args.gmapDir, args.threads)
    
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
        print("# An existing pickle file will be loaded to resume program operation.")
    
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        raise FileExistsError(f"'{pickleFile}' already exists within '{runDir}', but is not a file; " +
                              "move/delete/fix this to allow BINge to continue.")
    
    # ... or begin pre-external clustering BINge
    else:
        # Set up a bin collection structure for each genome
        collectionList = generate_bin_collections(locations.genomesDir, args.threads, args.isMicrobial)
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
    bingeResults.write_binge_clusters(outputFileName, clusterTypes="all")
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
    print("## File inputs:")
    print("# Genome targets:")
    genomeLinks = [ os.path.join(locations.genomesDir, f) for f in os.listdir(locations.genomesDir) if f.endswith(".fasta") ]
    annotLinks = []
    for genomeLink in genomeLinks:
        suffixNum = genomeLink.rsplit("genome", maxsplit=1)[1].split(".fasta")[0]
        annotLink = os.path.join(locations.genomesDir, f"genome{suffixNum}.gff3")
        if os.path.exists(annotLink):
            annotLinks.append(annotLink)
        else:
            annotLinks.append(None)
    genomeOrigins = [ str(Path(g).resolve()) for g in genomeLinks ]
    annotOrigins = [ str(Path(a).resolve()) if a is not None else None for a in annotLinks ]
    for gLink, gOrigin, aLink, aOrigin in zip(genomeLinks, genomeOrigins, annotLinks, annotOrigins):
        print(f"{gLink} -> {gOrigin}")
        if aOrigin is not None:
            print(f"    pre-seeded by: {aLink} -> {aOrigin}")
    
    print()
    
    print("# GFF3 annotations:")
    gff3Links = [ os.path.join(locations.gff3Dir, f) for f in os.listdir(locations.gff3Dir) if f.endswith(".gff3") ]
    gff3Origins = [ str(Path(g).resolve()) for g in gff3Links ]
    fastaLinks = [ os.path.join(locations.gff3Dir, f) for f in os.listdir(locations.gff3Dir) if f.endswith(".fasta") ]
    fastaOrigins = [ str(Path(f).resolve()) for f in fastaLinks ]
    if len(gff3Links) == 0:
        print("None")
    else:
        for gLink, gOrigin, fLink, fOrigin in zip(gff3Links, gff3Origins, fastaLinks, fastaOrigins):
            print(f"{gLink} -> {gOrigin}")
            print(f"    associated FASTA: {fLink} -> {fOrigin}")
    print()
    
    print("# Transcript FASTAs:")
    mrnaLinks = [ os.path.join(locations.txDir, f) for f in os.listdir(locations.txDir) if f.endswith(".mrna") ]
    otherSeqLinks = []
    for mrnaLink in mrnaLinks:
        suffixNum = mrnaLink.rsplit("transcriptome", maxsplit=1)[1].split(".mrna")[0]
        cdsLink = os.path.join(locations.txDir, f"transcriptome{suffixNum}.cds")
        protLink = os.path.join(locations.txDir, f"transcriptome{suffixNum}.aa")
        if os.path.exists(cdsLink) and os.path.exists(protLink):
            otherSeqLinks.append([cdsLink, protLink])
        else:
            otherSeqLinks.append(None)
    mrnaOrigins = [ str(Path(t).resolve()) for t in mrnaLinks ]
    otherSeqOrigins = [ (str(Path(value[0]).resolve()), str(Path(value[1]).resolve())) if value is not None else None for value in otherSeqLinks ]
    
    if len(mrnaLinks) == 0:
        print("None")
    else:
        ongoingCount = 0
        for mrnaLink, mrnaOrigin, otherSeqOrigin in zip(mrnaLinks, mrnaOrigins, otherSeqOrigins):
            if otherSeqOrigin is not None:
                cdsOrigin, protOrigin = otherSeqOrigin
                print(f"mRNA: {mrnaLink} -> {mrnaOrigin}")
                print(f"CDS: {cdsOrigin} -> {cdsOrigin}")
                print(f"AA: {protOrigin} -> {protOrigin}")
            else:
                suffixNum = mrnaLink.rsplit("transcriptome", maxsplit=1)[1].split(".mrna")[0]
                cdsFile = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.cds")
                aaFile = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.aa")
                print(f"mRNA: {mrnaLink} -> {mrnaOrigin}")
                print(f"CDS (predicted): {cdsFile}")
                print(f"AA (predicted): {aaFile}")
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
