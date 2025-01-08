#! python3
# BINge_post.py
# BIN Genes for Expression analyses - post-processing modules

# A collection of functions for downstream processing of BINge clustering,
# allowing you to 1) BLAST sequences against a database, 2) quantify reads
# using salmon, 3) filter clusters based on these results,
# and 4) extract representative sequences for each cluster also based on
# these results.

import os, argparse, sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts.Function_packages import ZS_BlastIO, ZS_MapIO

from modules.fasta_handling import ZS_SeqIO, FastaCollection
from modules.validation import validate_blast_args, validate_salmon_args, validate_filter_args, \
    validate_representatives_args, validate_dge_args, touch_ok
from modules.parsing import parse_equivalence_classes, parse_quants, parse_binge_clusters, \
    locate_read_files

# Define functions
def validate_args(args):
    # Derive the locations of the FASTA files
    args.fastaFiles = [
        os.path.join(args.fastaDir, f)
        for f in os.listdir(args.fastaDir)
        if f.endswith(".nucl")
    ]
    if len(args.fastaFiles) == 0:
        print(f"The -f directory does not contain any .nucl files.")
        print("Ideally, the -f argument is the location where BINge was run.")
        print("Please point me to a location containing .nucl files then try again.")
        quit()
    
    # Validate BLAST file if relevant
    if args.blastFile != None:
        if not os.path.isfile(args.blastFile):
            print(f'I am unable to locate the BLAST outfmt6 file ({args.blastFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate GFF3 / text file if relevant
    if args.annotationFile != None:
        if not os.path.isfile(args.annotationFile):
            print(f'I am unable to locate the annotation GFF3/text file ({args.annotationFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        # Check if we've received a text or GFF3 file
        with open(args.annotationFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                else:
                    l = line.rstrip("\r\n ")
                    sl = l.split("\t")
                    
                    if len(sl) == 9:
                        args.annotFileFormat = "gff3"
                    elif len(sl) == 1 and " " not in sl[0]:
                        args.annotFileFormat = "text"
                    else:
                        print(f"I was not able to determine the input file format of '{args.annotationFile}'")
                        print(f"Specifically, look at the first non-comment line '{l}'")
                        print("I expect it to have 9 tab-separated columns if a GFF3")
                        print("Or, I expect it to have no tab separation and no white space if it's a text file")
                        print("Make sure you have the right file in the right format, and try again.")
                        quit()

def parse_gff3_ids(gff3File):
    '''
    Simple function to iterate through a GFF3 file and hold onto any parent and
    subfeature IDs. This process should encompass genes and mRNAs as well as any
    other subfeature e.g., lnc_RNA, whilst skipping over exons and CDS'.
    '''
    # Setup for slim parsing of attributes
    def _format_attributes(attributes):
        "Code borrowed from ZS_GFF3IO"
        SLIM_ATTRIBUTES = ["id", "parent"]
        
        splitAttributes = []
        for a in attributes.split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        attributesDict = {
            splitAttributes[i]: splitAttributes[i+1]
            for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)
            if splitAttributes[i].lower() in SLIM_ATTRIBUTES
        }
        return attributesDict
    
    parentIDs = set()
    subIDs = set()
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Parse content lines
            if not line.startswith("#") and len(sl) == 9:
                featureType, attributes = sl[2], sl[8]
                attributesDict = _format_attributes(attributes)
                if featureType == "gene":
                    parentIDs.add(attributesDict["ID"])
                elif "Parent" in attributesDict and attributesDict["Parent"] in parentIDs:
                    subIDs.add(attributesDict["ID"])
    annotIDs = parentIDs.union(subIDs)
    return annotIDs

def parse_text_ids(textFile):
    '''
    Simply parses a text file for its sequence IDs. The input file is expected to be
    a newline-delimited list of IDs with no header.
    '''
    annotIDs = set()
    with open(textFile, "r") as fileIn:
        for line in fileIn:
            annotIDs.add(line.rstrip("\r\n "))
    return annotIDs

def parse_fasta_ids(fastaFiles):
    '''
    Simply parses one or more FASTA files for their sequence IDs. Returns a set
    of said IDs.
    '''
    fastaIDs = []
    for f in fastaFiles:
        with open(f, "r") as fileIn:
            for line in fileIn:
                if line.startswith(">"):
                    fastaIDs.append(line.rstrip("\r\n ").split(" ")[0])
    return set(fastaIDs)

def get_counts_cutoff_by_percentiles(binnedClusterDict, unbinnedClusterDict,
                                     salmonCollection, salmonFileFormat):
    '''
    Parses the sequence IDs out of the two cluster dictionaries, alongside a salmon
    EquivalenceClassCollection or QuantCollection (with format indicated), and calculates
    the nth percentile of read count values in the dataset.
    
    Parameters:
        binnedClusterDict -- a dictionary with format like:
                             {
                                 1: [ seqid1, seqid2, ... ],
                                 2: [...],
                                 ...
                             }
        unbinnedClusterDict -- a dictionary with the same format as binnedClusterDict.
        salmonCollection -- an EquivalenceClassCollection or QuantCollection.
        salmonFileFormat -- a string == "ec" if the collection is an EquivalenceClassCollection,
                            and "quant" if a QuantCollection
    '''
    # Get counts for each bin dictionary
    binnedCounts = [
        np.mean(salmonCollection.get_transcript_count(seqID))
        
        for seqIDs in binnedClusterDict.values()
        for seqID in seqIDs
        if (salmonFileFormat == "ec" and seqID in salmonCollection.ids)
        or (salmonFileFormat == "quant" and seqID in salmonCollection.quant)
    ]
    
    unbinnedCounts = [
        np.mean(salmonCollection.get_transcript_count(seqID))
        
        for seqIDs in unbinnedClusterDict.values()
        for seqID in seqIDs
        if (salmonFileFormat == "ec" and seqID in salmonCollection.ids)
        or (salmonFileFormat == "quant" and seqID in salmonCollection.quant)
    ]
    
    # If we have binnedCounts, derive a good cutoff from that
    percentiles = [20, 30, 40, 50, 60, 70, 80, 90, 99]
    if binnedCounts != []:
        # Find out where we first start seeing actual read alignments for binned dict
        for percentile in percentiles:
            binnedPercentileMean = np.percentile(binnedCounts, percentile)
            if binnedPercentileMean > 1: # == 1 read count on average
                break
        
        # Choose a percentile in the unbinned dict that is >= this number of counts
        for percentile in percentiles:
            unbinnedPercentileMean = np.percentile(unbinnedCounts, percentile)
            if unbinnedPercentileMean >= binnedPercentileMean:
                break
        
        # This becomes our countCutoff value
        "If not even the 99th percentile meets the cutoff, we just take it anyway"
        return unbinnedPercentileMean
    
    # If we do not have binnedCounts, just derive a cutoff de novo
    else:
        "This won't work as well, but is only a problem CD-HIT faces"
        for percentile in percentiles:
            unbinnedPercentileMean = np.percentile(unbinnedCounts, percentile)
            if unbinnedPercentileMean >= 1: # == 1 read count on average
                break
        
        # This becomes our countCutoff value (if it is >= 1)
        "If it doesn't even == 1, we just set 1 so this filter will NEVER save a cluster"
        return unbinnedPercentileMean if unbinnedPercentileMean >= 1 else 1

def merge_cluster_dicts(binnedDict, unbinnedDict, toDrop):
    '''
    Merges binnedDict and unbinnedDict together, after dropping any clusters
    noted in the toDrop set.
    
    Parameters:
        binnedDict -- a dictionary containing clusters with format like:
                          {
                              1: [ seqid1, seqid2, ... ],
                              2: [ ... ],
                              ...
                          }
        toFilterDict -- a dictionary with the same format as unfilteredDict.
        toDrop -- a set listing cluster numbers to exclude.
    Returns:
        clusterDict -- a new dictionary with clusters numbered in toDrop removed,
                       with a format like:
                       {
                           1: [ [seqid1, seqid2, ...], "binned" ],
                           2: [ [...], "unbinned"],
                           3: ...,
                           ...
                       }
    '''
    clusterDict = {}
    numClusters = 0
    for i in range(2):
        toFilterDict = [binnedDict, unbinnedDict][i]
        binType = "binned" if i == 0 else "unbinned"
        
        for clusterNum, seqIDs in toFilterDict.items():
            if not clusterNum in toDrop:
                numClusters += 1
                clusterDict[numClusters] = [seqIDs, binType]
    
    return clusterDict

def determine_if_1x_filter(transcriptCounts, seqIDs, transcriptRecords, readLength, sequenceType):
    '''
    Calculates if a cluster should be filtered on the basis of whether it has 1x
    coverage or not across all samples combined.
    
    Parameters:
        transcriptCounts -- a list containing sublists for each transcript in the
                            cluster, where sublists contain integer values indicating the
                            amount of read alignments made to that transcript from each
                            sample.
        seqIDs -- a list containing strings of transcript identifiers.
        transcriptRecords -- a pyfaidx.Fasta or FastaCollection which can be indexed to
                             retrieve all transcripts identified in transcriptIDs.
        readLength -- an integer indicating the average read length from sequencing.
        sequenceType -- a string in the list ["nucleotide", "protein"] indicating what type
                        of sequences we are working with.
    Returns:
        has1x -- a boolean indicating whether this cluster contains a transcript which meets
                 1x coverage criteria or not.
    '''
    # Get sequence lengths per transcript as nucleotides
    transcriptLengths = [ len(str(transcriptRecords[sid])) for sid in seqIDs ]
    transcriptLengths = [ tl * 3 if sequenceType == "protein" else tl for tl in transcriptLengths ]
    
    # Figure out how many reads would be required for 1x coverage
    readsNeededFor1x = [ tl / readLength for tl in transcriptLengths ]
    
    # Calculate if we've reached this threshold across samples
    has1x = any(
    [
        sum(tc) >= readsNeeded
        for tc, readsNeeded in zip(transcriptCounts, readsNeededFor1x)
    ])
    
    return has1x

def concatenate_sequences(sequenceFiles, outputFileName):
    '''
    Helper function to concatenate all sequence files into a single file, as needed
    when running BLAST or salmon read quantification.
    
    Parameters:
        sequenceFiles -- a list of strings indicating the locations of FASTA files to
                         concatenate.
        outputFileName -- a string indicating the location of the output file.
    '''
    if not os.path.exists(outputFileName) or not os.path.exists(outputFileName + ".ok"):
        print("# Concatenating all sequence files into a single file...")
        with open(outputFileName, "w") as fileOut:
            for fastaFile in sequenceFiles:
                with open(fastaFile, "r") as fileIn:
                    for line in fileIn:
                        fileOut.write(line)
                if not line.endswith("\n"):
                    fileOut.write("\n")
        touch_ok(outputFileName)

## Main
def main():
    blastDescription = """'blast' aims to query the sequences used during BINge clustering
    to a target database using MMseqs2. This information can be used to filter clusters or
    pick representatives of clusters."""
    
    salmonDescription = """'salmon' aims to quantify sequences used during BINge clustering
    with reads from one or more FASTQ files. This information can be used to filter clusters
    or pick representatives of clusters. It also allows for their downstream use in DGE
    analysis."""
    
    filterDescription = """'filter' aims to remove low-quality clusters from a BINge
    clustering analysis on the basis of various biological evidence. These can include
    sequence length, BLAST evidence, read alignment, and reference genome annotations.
    Unbinned clusters (or all clusters if providing --alsoFilterBinned) will be retained
    if they contain sequence(s) which meet any of the below criteria:
    1) sequence comes from a reference GFF3, or
    2) has a significant BLAST E-value, or
    3) has read alignment > cut-off value from more than 1 sample (if multiple
    Salmon files are given), or
    4) has an ORF >= cut-off length.
    
    As an example, if filter 2 passes i.e., a sequence in a cluster has a significant
    BLAST hit, that cluster will NOT be filtered even if its ORF length is less than
    that cut-off.
    
    The output is modified BINge cluster TSV file with failing clusters removed. Cluster
    numbers will be consistent with the original unfiltered file.
    """
    
    representativesDescription = """'representatives' aims to pick a single representative
    sequence from each cluster in a BINge clustering analysis. This is useful for downstream
    analyses where you want to annotate genes or perform differential gene expression.
    Representatives are picked based on biological evidence in order of importance:
    (is a reference organism sequence) > (BLAST E-value) > (number of reads aligned) >
    (sequence length). If you do not indicate the use of BLAST or salmon read alignment
    evidence, this script will pick a  representative on the basis of sequence length
    (longest being best). The output is a FASTA file containing these representatives.
    The ID for each sequence follows a format like '>Cluster-1 representative=transcript_92'.
    """
    
    dgeDescription = """TBD"""
    
    mainDescription = """%(prog)s provides the ability to 'filter' clusters from a BINge
    clustering analysis, as well as pick 'representatives' for each cluster. Data required
    to do this (i.e., BLAST and salmon read quantification) is facilitated by the 'blast' and
    'salmon' options, respectively. Note that the term 'BLAST' is being used informally;
    the MMseqs2 program is used for any sequence similarity searches."""
    
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-d", dest="workingDirectory",
                   required=True,
                   help="Specify the location where clustering has been performed")
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=mainDescription)
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    bparser = subparsers.add_parser("blast",
                                    parents=[p],
                                    add_help=False,
                                    help="Run BLAST/MMseqs2 query against a database",
                                    description=blastDescription)
    bparser.set_defaults(func=bmain)
    
    sparser = subparsers.add_parser("salmon",
                                    parents=[p],
                                    add_help=False,
                                    help="Run salmon read quantification against sequences",
                                    description=salmonDescription)
    sparser.set_defaults(func=smain)
    
    fparser = subparsers.add_parser("filter",
                                    parents=[p],
                                    add_help=False,
                                    help="Run filtering to remove low-quality clusters",
                                    description=filterDescription)
    fparser.set_defaults(func=fmain)
    
    rparser = subparsers.add_parser("representatives",
                                    parents=[p],
                                    add_help=False,
                                    help="Pick representative sequences for each cluster",
                                    description=representativesDescription)
    rparser.set_defaults(func=rmain)
    
    dparser = subparsers.add_parser("dge",
                                    parents=[p],
                                    add_help=False,
                                    help="Format inputs for downstream DGE analysis",
                                    description=dgeDescription)
    dparser.set_defaults(func=dmain)
    
    # BLAST-subparser arguments
    bparser.add_argument("-t", dest="targetFile",
                         required=True,
                         help="""Specify the location of a FASTA file to query against.""")
    bparser.add_argument("-s", dest="sequenceType",
                         required=True,
                         choices=["protein", "nucleotide"],
                         help="""Specify whether your database contains 'protein' or 'nucleotide'
                         sequences.""")
    ## Optional (program locations)
    bparser.add_argument("--mms2", dest="mms2Exe",
                         required=False,
                         help="""Specify the location of the MMseqs2 executable if not locateable
                         in your PATH variable.""",
                         default=None)
    ## Optional (parameters)
    bparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run search with
                         (default==1)""",
                         default=1)
    
    # Salmon-subparser arguments
    sparser.add_argument("-r", dest="readsDir",
                         required=True,
                         nargs="+",
                         help="""Specify one or more locations containing FASTQ files
                         for read quantification.""")
    sparser.add_argument("-s", dest="readsSuffix",
                         required=True,
                         help="""Indicate the suffix of your FASTQ files (e.g., '.fastq',
                         '.fq', '.fq.gz'); if your reads are paired, this suffix must
                         immediately follow the read number (e.g., '_1.fastq', '_2.fastq');
                         all files with this suffix in locations provided
                         to -r will be used.""")
    ## Optional (flags)
    sparser.add_argument("--singleEnd", dest="singleEnd",
                         required=False,
                         action="store_true",
                         help="""Set this flag if your reads are single-ended; if they are not,
                         read pairs will be automatically determined by shared file prefix.""",
                         default=False)
    ## Optional (program locations)
    sparser.add_argument("--salmon", dest="salmonExe",
                         required=False,
                         help="""Specify the location of the salmon executable if not locateable
                         in your PATH variable.""",
                         default=None)
    ## Optional (parameters)
    sparser.add_argument("--threads", dest="threads",
                         required=False,
                         type=int,
                         help="""Optionally, specify how many threads to run search with
                         (default==1)""",
                         default=1)
    
    # Filter-subparser arguments
    fparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis folder to view by its hash; if not provided,
                         the most recent analysis folder will be viewed""",
                         default="most_recent")
    ## Optional (flags)
    fparser.add_argument("--justDropUnbinned", dest="justDropUnbinned",
                         required=False,
                         action="store_true",
                         help="""Set this flag to remove all unbinned clusters but keep all binned
                         clusters; this supercedes all other filters.""",
                         default=False)
    fparser.add_argument("--alsoFilterBinned", dest="filterBinned",
                         required=False,
                         action="store_true",
                         help="""Optionally, extend filtering to clusters that were binned. By default,
                         binned clusters are assumed to be good since they've already met sequence
                         similarity criteria to the reference genome hence filtering is applied only
                         to unbinned clusters.""",
                         default=False)
    fparser.add_argument("--useGFF3", dest="useGFF3",
                         required=False,
                         action="store_true",
                         help="""Set this flag to keep clusters if they contain any
                         sequences found within provided GFF3 annotations (-ig during 'initialise').""",
                         default=False)
    fparser.add_argument("--useBLAST", dest="useBLAST",
                         required=False,
                         action="store_true",
                         help="""Set this flag to keep clusters if they contain any sequences with
                         significant BLAST hits.""",
                         default=False)
    fparser.add_argument("--useSalmon", dest="useSalmon",
                         required=False,
                         action="store_true",
                         help="""Set this flag to keep clusters if they contain any sequences with
                         read alignments above 1x coverage threshold.""",
                         default=False)
    ## Optional (parameters)
    fparser.add_argument("--evalue", dest="evalue",
                         required=False,
                         type=float,
                         help="""If using --keepBLAST, specify the E-value cut-off for BLAST hits
                         to be considered significant (default==1e-10).""",
                         default=1e-10)
    fparser.add_argument("--readLength", dest="readLength",
                         required=False,
                         type=int,
                         help="""If specifying --keepSalmon, indicate the length of your reads in
                         basepairs (default==150; just consider the single end length)""",
                         default=150)
    fparser.add_argument("--length", dest="minimumLength",
                         required=False,
                         type=int,
                         help="""Specify the minimum length of an ORF (in amino acids) that must
                         be met for a cluster to pass filtration (if it hasn't passed any other
                         filter checks already); default==0""",
                         default=0)
    
    # Representatives-subparser arguments
    rparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis folder to view by its hash; if not provided,
                         the most recent analysis folder will be viewed""",
                         default="most_recent")
    ## Optional (flags)
    rparser.add_argument("--useGFF3", dest="useGFF3",
                         required=False,
                         action="store_true",
                         help="""Set this flag to pick a representative if it comes from a GFF3
                         annotation file (-ig during 'initialise').""",
                         default=False)
    rparser.add_argument("--useBLAST", dest="useBLAST",
                         required=False,
                         action="store_true",
                         help="""Set this flag to pick a representative that has a
                         significant BLAST hit.""",
                         default=False)
    rparser.add_argument("--useSalmon", dest="useSalmon",
                         required=False,
                         action="store_true",
                         help="""Set this flag to pick a representative that has the most
                         read alignments.""",
                         default=False)
    ## Optional (parameters)
    rparser.add_argument("--evalue", dest="evalue",
                         required=False,
                         type=float,
                         help="""If using --keepBLAST, specify the E-value cut-off for BLAST hits
                         to be considered significant (default==1e-10).""",
                         default=1e-10)
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "blast":
        print("## BINge_post.py - MMseqs2 query ##")
        validate_blast_args(args) # sets args.sequenceFiles
        bmain(args)
    elif args.mode == "salmon":
        print("## BINge_post.py - Salmon read quantification ##")
        validate_salmon_args(args)
        smain(args)
    elif args.mode == "filter":
        print("## BINge_post.py - Filter clusters ##")
        validate_filter_args(args) # sets args.runDirName, args.bingeFile
        fmain(args)
    elif args.mode == "representatives":
        print("## BINge_post.py - Representatives selection ##")
        validate_representatives_args(args) # sets args.runDirName, args.bingeFile
        rmain(args) # sets args.bingeFile
    elif args.mode == "dge":
        print("## BINge_post.py - DGE preparation ##")
        validate_dge_args(args) # sets args.runDirName, args.bingeFile
        dmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args):
    # Set up BLAST directory
    blastDir = os.path.join(args.workingDirectory, "blast")
    os.makedirs(blastDir, exist_ok=True)
    
    tmpDir = os.path.join(blastDir, "tmp")
    os.makedirs(tmpDir, exist_ok=True)
    
    # Figure out if we've already run BLAST and exit if so
    outputFileName = os.path.join(blastDir, f"MMseqs2_results.tsv")
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"A MMseqs2 output file already exists within '{blastDir}'; " +
                              "this program will not overwrite an existing file. " +
                              "To resume or overwrite, move/rename/delete the file " + 
                              f"'{outputFileName}' then try again.")
    
    # Concatenate all FASTA files into a single file
    sequenceSuffix = ".cds" if args.sequenceType == "nucleotide" else ".aa"
    concatFileName = os.path.join(blastDir, f"concatenated{sequenceSuffix}")
    concatenate_sequences(args.sequenceFiles, concatFileName)
    
    # Establish the MMseqs2 query and target databases
    queryDB = ZS_BlastIO.MM_DB(concatFileName, os.path.dirname(args.mms2Exe), tmpDir,
                               args.sequenceType, args.threads)
    queryDB.generate()
    queryDB.index()
    
    targetDB = ZS_BlastIO.MM_DB(args.targetFile, os.path.dirname(args.mms2Exe), tmpDir,
                                args.sequenceType, args.threads)
    targetDB.generate()
    targetDB.index()
    
    # Run MMseqs2 search
    mms2 = ZS_BlastIO.MMseqs(queryDB, targetDB, os.path.dirname(args.mms2Exe), tmpDir)
    mms2.threads = args.threads
    mms2.evalue = 1 # weak threshold, leave it to 'filter' to decide on a good one
    mms2.mmseqs(outputFileName, force=True) # force=True to overwrite any existing file if '.ok' is missing
    touch_ok(outputFileName)
    
    print("BLAST search complete!")

def smain(args):
    # Set up salmon directory
    salmonDir = os.path.join(args.workingDirectory, "salmon")
    os.makedirs(salmonDir, exist_ok=True)
    
    # Locate read files
    forwardReads, reverseReads, sampleNames = locate_read_files(args.readsDir, args.readsSuffix,
                                                               args.singleEnd)
    
    # Concatenate all FASTA files into a single file
    concatFileName = os.path.join(salmonDir, f"concatenated.cds") # map to CDS
    concatenate_sequences(args.sequenceFiles, concatFileName)
    
    # Establish the salmon target database
    targetDB = ZS_MapIO.Salmon_DB(concatFileName, os.path.dirname(args.salmonExe),
                                  args.threads)
    if not targetDB.index_exists():
        targetDB.index()
    
    # Run salmon read quantification
    for i, sampleName in enumerate(sampleNames):
        # Establish the salmon object
        if args.singleEnd:
            salmon = ZS_MapIO.Salmon([forwardReads[i]], targetDB, os.path.dirname(args.salmonExe),
                                     args.threads)
        else:
            salmon = ZS_MapIO.Salmon([forwardReads[i], reverseReads[i]], targetDB, os.path.dirname(args.salmonExe),
                                     args.threads)
        
        # Set up sample directory
        sampleDir = os.path.join(salmonDir, sampleName)
        os.makedirs(sampleDir, exist_ok=True)
        
        # Run salmon (if not already run)
        if not os.path.exists(sampleDir + ".ok"):
            print(f"# Running salmon quantification for '{sampleName}'...")
            salmon.quant(sampleDir)
            touch_ok(sampleDir + ".ok")
    
    print("Salmon quant complete!")

def fmain(args):
    # Set up filter directory for this run
    filterDir = os.path.join(args.workingDirectory, "filtered")
    os.makedirs(filterDir, exist_ok=True)
    
    runDir = os.path.join(filterDir, args.runDirName)
    os.makedirs(runDir, exist_ok=True)
    
    mostRecentDir = os.path.join(filterDir, "most_recent")
    if os.path.exists(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(runDir, mostRecentDir)
    
    # Figure out if we've already filtered this run and exit if so
    outputFileName = os.path.join(runDir, f"BINge_clustering_result.filtered.tsv")
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"A filtered BINge output file already exists within '{runDir}'; " +
                              "this program will not overwrite an existing file. " +
                              "To resume or overwrite, move/rename/delete this file then try again.")
    
    # Parse the BINge cluster file
    binnedDict = parse_binge_clusters(args.bingeFile, "binned")
    unbinnedDict = parse_binge_clusters(args.bingeFile, "unbinned")
    
    # If we're just dropping unbinned clusters, we can do that now
    if args.justDropUnbinned:
        with open(outputFileName, "w") as fileOut:
            # Write header
            fileOut.write("#BINge clustering information file\n")
            fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
            
            # Write content lines
            for clusterNum, seqIDs in binnedDict.items():
                for seqID in seqIDs:
                    fileOut.write(f"{clusterNum}\t{seqID}\tbinned\n")
        return # exit the program here
    
    # No testing done beyond this point yet
    return
    
    # Load transcripts into memory for quick access
    transcriptRecords = FastaCollection(args.fastaFiles)
    
    # Parse BLAST results (if relevant)
    if args.blastFile != None:
        blastResults = ZS_BlastIO.BLAST_Results(args.blastFile)
        blastResults.evalue = args.evalue
        blastResults.num_hits = 1 # only need to keep the best hit for each sequence
        blastResults.parse()
        blastDict = blastResults.results
    else:
        blastDict = {}
    
    # Parse annotation file (if relevant)
    if args.annotationFile != None:
        if args.annotFileFormat == "gff3":
            annotIDs = parse_gff3_ids(args.annotationFile)
        else:
            annotIDs = parse_text_ids(args.annotationFile)
    else:
        annotIDs = set()
    
    # Parse filter FASTA files (if relevant)
    if args.filterFiles != []:
        filterIDs = parse_fasta_ids(args.filterFiles)
    else:
        filterIDs = None
    
    # Parse salmon counts (if relevant)
    if args.salmonFiles != []:
        sampleNames = [ f"{i}" for i in range(len(args.salmonFiles))] # sample names don't matter
        
        if args.salmonFileFormat == "ec":
            salmonCollection = parse_equivalence_classes(args.salmonFiles, sampleNames)
        else:
            salmonCollection = parse_quants(args.salmonFiles, sampleNames)
    else:
        salmonCollection = None
    
    # Get an appropriate read count value to use as a cut-off
    """Note: we want this to be somewhat strict since each cut-off SAVES a sequence,
    so if this is the only thing that passes on a cluster, we want to be confident
    that it's actually good."""
    countCutoff = get_counts_cutoff_by_percentiles(
        binnedDict, unbinnedDict,
        salmonCollection, args.salmonFileFormat
    )
    
    # Perform filtration of clusters with available evidence
    toDrop = set()
    for index, toFilterDict in enumerate([binnedDict, unbinnedDict]):
        for clusterID, seqIDs in toFilterDict.items():
            # Check 0: FILTER if it does not contain a filter ID
            if args.filterFiles != []:
                if not any([ seqID in filterIDs for seqID in seqIDs ]):
                    toDrop.add(clusterID)
                    continue # since we're filtering, we toDrop it then continue
            
            # Skip now if we don't want to filter binned clusters
            if index == 0 and (not args.filterBinned):
                continue
            
            # Check 1: Retain if it has a reference sequence
            if any([ seqID in annotIDs for seqID in seqIDs ]):
                continue
            
            # Handle read alignment values
            if salmonCollection != None:
                transcriptCounts = [
                    salmonCollection.get_transcript_count(seqID)
                    for seqID in seqIDs
                    if (args.salmonFileFormat == "ec" and seqID in salmonCollection.ids)
                    or (args.salmonFileFormat == "quant" and seqID in salmonCollection.quant)
                ]
                
                # Check 2: FILTER if it lacks 1x coverage
                if args.require1x:
                    has1x = determine_if_1x_filter(transcriptCounts, seqIDs,
                                   transcriptRecords, args.readLength,
                                   "protein" if args.sequenceFormat == "protein"
                                   else "nucleotide")
                    
                    if not has1x:
                        toDrop.add(clusterID)
                        continue # since we're filtering, we toDrop it then continue
                
                # Check 3: Retain if any have good read alignment in >1 sample
                if len(salmonCollection.samples) > 1:
                    if any(
                    [
                        sum([ countValue > countCutoff for countValue in tc ]) > 1
                        for tc in transcriptCounts
                    ]):
                        continue
                
                # Check 3 (alt): Retain if it has good read alignment in the 1 sample
                else:
                    if any(
                    [
                        tc[0] > countCutoff
                        for tc in transcriptCounts
                    ]):
                        continue
            
            # Check 4: Retain if any have a significant E-value
            if any([ seqID in blastDict for seqID in seqIDs ]): # if in blastDict, it has a good E-value
                continue
            
            # Check 5: Retain if the sequence is long enough
            "Note: we include stop codons in the length calculations here"
            if args.sequenceFormat == "protein":
                if any(
                [
                    len(str(transcriptRecords[seqID])) >= args.minimumLength
                    for seqID in seqIDs
                ]):
                    continue
            elif args.sequenceFormat == "cds":
                if any(
                [
                    ( len(str(transcriptRecords[seqID])) / 3 ) >= args.minimumLength
                    for seqID in seqIDs
                ]):
                    continue
            else:
                seqs = [ ZS_SeqIO.FastASeq("x", str(transcriptRecords[seqID])).get_translation(True)[0] for seqID in seqIDs ]
                if any(
                [
                    len(seq[seq.find("M"):]) >= args.minimumLength
                    for seq in seqs
                ]):
                    continue
            
            # If we make it here, this cluster should be filtered
            toDrop.add(clusterID)
    
    # Merge cluster dicts back together
    clusterDict = merge_cluster_dicts(binnedDict, unbinnedDict, toDrop)
    
    # Write filtered clusters to file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines
        for clusterNum, value in clusterDict.items():
            clusterIDs, binType = value
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum}\t{seqID}\t{binType}\n")
    
    # Print some statistics for the user
    print("# BINge_post.py 'filter' statistics:")
    print(f" > used 1x criteria to filter clusters = {args.require1x}")
    print(f" > number of reads = {countCutoff} required for read count rescue")
    print(f" > removed {len(toDrop)} clusters")
    
    print("Filtering complete!")

def rmain(args):
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeFile,
        "all" if args.onlyBinned == False else "binned"
    )
    
    # No testing done beyond this point yet
    return
    
    # Load transcripts into memory for quick access
    transcriptRecords = FastaCollection(args.fastaFiles) # args.fastaFiles set by validate_args()
    
    # Parse BLAST results (if relevant)
    if args.blastFile != None:
        blastDict = parse_blast_outfmt6(args.blastFile, args.evalue)
    else:
        blastDict = {}
    
    # Parse annotation file (if relevant)
    if args.annotationFile != None:
        if args.annotFileFormat == "gff3":
            annotIDs = parse_gff3_ids(args.annotationFile)
        else:
            annotIDs = parse_text_ids(args.annotationFile)
    else:
        annotIDs = set()
    
    # Parse salmon counts (if relevant)
    if args.salmonFiles != []:
        sampleNames = [ f"{i}" for i in range(len(args.salmonFiles))] # sample names don't matter
        
        if args.salmonFileFormat == "ec":
            salmonCollection = parse_equivalence_classes(args.salmonFiles, sampleNames)
        else:
            salmonCollection = parse_quants(args.salmonFiles, sampleNames)
    else:
        salmonCollection = None
    
    # Write cluster representatives to file
    with open(args.outputFileName, "w") as fileOut:
        for clusterID, seqIDs in clusterDict.items():
            # Handle single sequence clusters (by just choosing the sequence)
            if len(seqIDs) == 1:
                seqID = seqIDs[0]
            
            # Handle multisequence clusters
            else:
                evidenceLists = []
                
                # Store all forms of evidence we have at hand for each sequence member
                numRepsInCluster = sum([ seqID in annotIDs for seqID in seqIDs ])
                for seqID in seqIDs:
                    # Tolerantly handle non-existing sequences
                    """This is only allowed because bin seeding via GFF3 may lead to binned sequences
                    which do not exist in the FASTA files. We can ignore these as their mRNA counterparts
                    should have been binned here as well."""
                    if not seqID in transcriptRecords:
                        if args.beTolerant:
                            continue
                        else:
                            print(f"Error: sequence ID '{seqID}' not found in any input FASTA files!")
                            print("--be_tolerant was not specified so the program will end now.")
                            quit()
                    
                    # Tolerantly handle counts that can be filtered by salmon
                    try:
                        counts = sum(salmonCollection.get_transcript_count(seqID))
                    except: # this happens if Salmon filtered something out
                        counts = 0
                    
                    # Now store the evidence
                    """We prioritise annotID values ONLY if there's a single one in the cluster;
                    if there are 2 or more, BINge may have merged fragmented gene clusters together
                    in which case selecting the original annotID sequences would end up ignoring the
                    'fix' that we've attempted to bring about. Moreover, if we've merged bins together
                    through self linking, it might be misleading to always present the original ID as
                    the representative since we've clustered more than one original and hence it becomes
                    a bit arbitrary to pick one above the other regardless of E-value BLAST score."""
                    if numRepsInCluster <= 1:
                        thisEvidenceList = [
                            1 if seqID in annotIDs else 0,
                            blastDict[seqID] if seqID in blastDict else 0, # stores bitscore
                            counts,
                            len(str(transcriptRecords[seqID])),
                            seqID
                        ]
                    else:
                        thisEvidenceList = [
                            0,
                            blastDict[seqID] if seqID in blastDict else 0, # stores bitscore
                            counts,
                            len(str(transcriptRecords[seqID])),
                            seqID
                        ]
                    evidenceLists.append(thisEvidenceList)
                
                # Raise error if we skipped all the IDs in this bin
                if args.beTolerant and len(evidenceLists) == 0:
                    print("--be_tolerant behaviour led to the discovery of an empty bin!")
                    print("Currently I will not handle this situation. You should create " + 
                          "extra FASTA(s) containing parent gene identifiers and use those as input.")
                    print("Program will exit now.")
                    quit()
                
                # Sort our evidence lists in a way where the first value is the best
                evidenceLists.sort(
                    key = lambda x: (
                        -x[0], -x[1], -x[2], -x[3]
                    )
                )
                
                # Choose the top sequence
                seqID = evidenceLists[0][-1]
            
            # Write to file
            seq = str(transcriptRecords[seqID])
            fastaString = format_representative(clusterID, seqID, seq)
            fileOut.write(fastaString)

    print("Representative picking complete!")

def dmain(args):
    pass

if __name__ == "__main__":
    main()
