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
from goatools import obo_parser

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts.Function_packages import ZS_BlastIO, ZS_MapIO

from modules.locations import Locations
from modules.fasta_handling import ZS_SeqIO, FastaCollection
from modules.salmon import SalmonQC
from modules.validation import validate_blast_args, validate_salmon_args, validate_filter_args, \
    validate_representatives_args, validate_dge_args, validate_annotate_args, touch_ok
from modules.parsing import parse_equivalence_classes, parse_quants, parse_binge_clusters, \
    parse_gff3_ids, locate_read_files, parse_binge_representatives

# Define functions
def get_counts_cutoff_by_percentiles(binnedClusterDict, unbinnedClusterDict,
                                     quantCollection):
    '''
    Parses the sequence IDs out of the two cluster dictionaries, alongside a salmon
    QuantCollection, and calculates the nth percentile of read count values in the
    dataset.
    
    Parameters:
        binnedClusterDict -- a dictionary with format like:
                             {
                                 1: [ seqid1, seqid2, ... ],
                                 2: [...],
                                 ...
                             }
        unbinnedClusterDict -- a dictionary with the same format as binnedClusterDict.
        quantCollection -- a QuantCollection.
    '''
    # Get counts for each bin dictionary
    binnedCounts = [
        np.mean(quantCollection.get_transcript_count(seqID))
        
        for seqIDs in binnedClusterDict.values()
        for seqID in seqIDs
        if seqID in quantCollection.quant
    ]
    
    unbinnedCounts = [
        np.mean(quantCollection.get_transcript_count(seqID))
        
        for seqIDs in unbinnedClusterDict.values()
        for seqID in seqIDs
        if seqID in quantCollection.quant
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
                           3: [ [...], "unbinned"],
                           12: ...,
                           ...
                       }
    '''
    clusterDict = {}
    for i in range(2):
        toFilterDict = [binnedDict, unbinnedDict][i]
        binType = "binned" if i == 0 else "unbinned"
        
        for clusterNum, seqIDs in toFilterDict.items():
            if not clusterNum in toDrop:
                clusterDict[clusterNum] = [seqIDs, binType]
    
    return clusterDict

def determine_if_1x_filter(transcriptCounts, seqIDs, transcriptRecords, readLength, sequenceType="nucleotide"):
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

def format_representative(clusterNum, representativeID, representativeSeq):
    '''
    Generates a string representation of a FASTA record with features to indicate
    the representative sequence.
    
    Parameters:
        clusterNum -- an int or string digit identifying the cluster.
        representativeID -- a string of the sequence ID of the representative of
                            this cluster.
        representativeSeq -- a string of the sequence itself for the representative
                             of this cluster.
    Returns:
        fastaString -- a string of the representative sequence formatted for writing
                       to file.
    '''
    return f">cluster-{clusterNum} representative={representativeID}\n{representativeSeq}\n"

def write_tx2gene(tx2geneFile, clusterDict):
    '''
    Generates a tx2gene file for summarisation of salmon counts to gene/cluster level
    from mapping done to individual transcripts/sequences.
    
    Parameters:
        tx2geneFile -- a string indicating the location of the output file.
        clusterDict -- a dictionary with format like:
                       {
                           1: [ seqid1, seqid2, ... ],
                           2: [ ... ],
                           ...
                       }
    '''
    with open(tx2geneFile, "w") as fileOut:
        fileOut.write("TXNAME\tGENEID\n")
        for clusterNum, seqIDs in clusterDict.items():
            for seqID in seqIDs:
                fileOut.write(f"{seqID}\tcluster-{clusterNum}\n")

def write_salmonQC(salmonDirs, salmonQCFile, clusterDict):
    '''
    Generates a salmonQC file for summarisation of salmon counts to gene/cluster level
    from mapping done to individual transcripts/sequences.
    
    Parameters:
        salmonDirs -- a list of strings indicating the locations of the salmon output directories.
        salmonQCFile -- a string indicating the location of the output file.
        clusterDict -- a dictionary with format like:
                       {
                           1: [ seqid1, seqid2, ... ],
                           2: [ ... ],
                           ...
                       }
    '''
    # Generate the SalmonQC object
    salmonQC = SalmonQC(salmonDirs)
    
    # Calculate adjusted mapping percentages
    adjustedQCDict = salmonQC.get_clustered_qc(clusterDict)
    
    # Write output
    with open(salmonQCFile, "w") as fileOut:
        fileOut.write("sample\traw_mapping_percent\tfiltered_mapping_percent\n")
        for sample, rawMapPct in salmonQC.mapPct.items():
            _, filteredMapPct = adjustedQCDict[sample] # drop the retainedPct
            fileOut.write(f"{sample}\t{rawMapPct}\t{filteredMapPct:.4f}\n")

def write_samples(salmonDirs, sampleFile):
    '''
    Generates a file which simply lists the samples used in the salmon quantification.
    
    Parameters:
        salmonDirs -- a list of strings indicating the locations of the salmon output directories.
        sampleFile -- a string indicating the location of the output file.
    '''
    # Get sample names from the salmon directories
    sampleNames = [ os.path.basename(salmonDir) for salmonDir in salmonDirs ]
    
    # Write output
    with open(sampleFile, "w") as fileOut:
        fileOut.write("sample\n")
        for sample in sampleNames:
            fileOut.write(f"{sample}\n")

def write_r_script(salmonFolder, sampleFile, tx2geneFile, rScriptFile):
    '''
    Prints ( / presents) instructions for how to load in the resulting files
    from this program into R.
    
    Parameters:
        salmonFolder -- a string indicating the location of the salmon output directories.
        sampleFile -- a string indicating the location of the samples file which lists
                       each sample name (which should match the salmon output directories).
        tx2geneFile -- a string indicating the location of the tx2gene file.
        rScriptFile -- a string indicating the location of the output R script.
    '''
    # Format text of R script
    rScriptText = f'''library(DESeq2)
library(tximport)
library(readr)

# Locate files
QUANT_FOLDER <- "{salmonFolder}"
SAMPLES_FILE <- "{sampleFile}"
TX2GENE_FILE <- "{tx2geneFile}"

# Parse sample names
samples <- read.table(file=SAMPLES_FILE, header = TRUE)[,1]

# Locate salmon quant files
quantFiles <- file.path(QUANT_FOLDER, samples, "quant.sf")
names(quantFiles) <- samples

# Parse tx2gene file
tx2gene <- read_delim(Tx2GENE_FILE, delim="\\t")

# tximport data
txi <- tximport(quantFiles, type="salmon", tx2gene=tx2gene)

# Load into DESeqDataSet object
# dds <- DESeqDataSetFromTximport(txi = txi, colData = < ... >, design = < ~ ... >)
'''
    # Write to file
    with open(rScriptFile, "w") as fileOut:
        fileOut.write(rScriptText)

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
    
    dgeDescription = """'dge' aims to streamline the preparation of data for downstream
    differential gene expression (DGE) analysis. You must have already performed clustering
    and 'salmon' herein to use this script, although you probably should also use 'filter'
    to remove low-quality or unbinned clusters. This script will generate an R script which
    will load the necessary data for DGE analysis with DESeq2 at your own discretion."""
    
    annotateDescription = """'annotate' aims to produce an annotation table using BLAST results
    for the representative sequences of each cluster. BLAST must have been run on the sequences
    with a UniRef database, and the 'idmapping_selected.tab' file from UniProtKB must be provided.
    The output is a TSV file indicating BLAST results as well as GO terms."""
    
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
                                    help="Run MMseqs2 query against a database",
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
    
    aparser = subparsers.add_parser("annotate",
                                    parents=[p],
                                    add_help=False,
                                    help="Produce an annotation table for clusters",
                                    description=annotateDescription)
    aparser.set_defaults(func=dmain)
    
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
                         help="""Specify the analysis folder to filter; if not provided,
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
                         sequences found within provided GFF3 annotations (-ig and/or
                         -t during 'initialise').""",
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
                         filter checks already); default==300 (bp)""",
                         default=300)
    
    # Representatives-subparser arguments
    rparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis or filter folder to obtain representatives from;
                         if not provided, the most recent analysis folder will be viewed""",
                         default="most_recent")
    ## Optional (flags)
    rparser.add_argument("--useGFF3", dest="useGFF3",
                         required=False,
                         action="store_true",
                         help="""Set this flag to pick a representative if it comes from a GFF3
                         annotation file (-ig and/or -t during 'initialise').""",
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
    
    # DGE-subparser arguments
    dparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis or filter folder to run DGE on;
                         if not provided, the most recent analysis folder will be viewed""",
                         default="most_recent")
    
    # Annotate-subparser arguments
    aparser.add_argument("-id", dest="idmappingFile",
                         required=True,
                         help="""Specify the location of the 'idmapping_selected.tab' file""")
    aparser.add_argument("-io", dest="oboFile",
                         required=True,
                         help="""Specify the location of a 'go.obo' file""")
    aparser.add_argument("--analysis", dest="analysisFolder",
                         required=False,
                         help="""Specify the analysis or filter folder to annotate clusters for;
                         if not provided, the most recent analysis folder will be viewed""",
                         default="most_recent")
    ## Optional (parameters)
    aparser.add_argument("--evalue", dest="evalue",
                         required=False,
                         type=float,
                         help="""Specify the E-value cut-off for BLAST hits
                         to be considered significant (default==1e-10).""",
                         default=1e-10)
    aparser.add_argument("--numhits", dest="numHits",
                         required=False,
                         type=int,
                         help="""Specify the number of hits to report for each sequence
                         (default==10).""",
                         default=10)
    aparser.add_argument("--largeTable", dest="largeTable",
                         required=False,
                         action="store_true",
                         help="""Set this flag to create a table with full outfmt6 fields
                         (e.g., gap opens, query start, etc.) rather than just identity,
                         E-value, and bitscore.""",
                         default=False)
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "blast":
        print("## BINge_post.py - MMseqs2 query ##")
        locations = validate_blast_args(args) # sets sequenceFiles
        bmain(args, locations)
    elif args.mode == "salmon":
        print("## BINge_post.py - Salmon read quantification ##")
        locations = validate_salmon_args(args) # sets sequenceFiles
        smain(args, locations)
    elif args.mode == "filter":
        print("## BINge_post.py - Filter clusters ##")
        locations = validate_filter_args(args) # sets sequenceFiles, runDirName, bingeFile, blastFile, gff3Files, salmonFiles
        fmain(args, locations)
    elif args.mode == "representatives":
        print("## BINge_post.py - Representatives selection ##")
        locations = validate_representatives_args(args) # sets sequenceFiles, runDirName, bingeFile, blastFile, gff3Files, salmonFiles
        rmain(args, locations) # sets args.bingeFile
    elif args.mode == "dge":
        print("## BINge_post.py - DGE preparation ##")
        locations = validate_dge_args(args) # sets runDirName, bingeFile, salmonFiles
        dmain(args, locations)
    elif args.mode == "annotate":
        print("## BINge_post.py - Cluster annotation ##")
        locations = validate_annotate_args(args) # sets runDirName, bingeFile, blastFile, targetFile, databaseTag
        amain(args, locations)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args, locations):
    # Set up BLAST directory
    os.makedirs(locations.blastDir, exist_ok=True)
    
    tmpDir = os.path.join(locations.blastDir, "tmp")
    os.makedirs(tmpDir, exist_ok=True)
    
    # Symlink to target file
    targetLink = os.path.join(locations.blastDir, locations.targetFile)
    if os.path.exists(targetLink) or os.path.islink(targetLink):
        os.unlink(targetLink)
    os.symlink(args.targetFile, targetLink)
    
    # Figure out if we've already run BLAST and exit if so
    outputFileName = os.path.join(locations.blastDir, locations.blastFile)
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"The MMseqs2 output file '{locations.blastFile}' already exists " +
                              f"within '{locations.blastDir}'; this program will not overwrite an existing " + 
                              "file. To resume program operation, move/rename/delete the file " + 
                              "then try again.")
    
    # Concatenate all FASTA files into a single file
    sequenceSuffix = ".cds" if args.sequenceType == "nucleotide" else ".aa"
    concatFileName = os.path.join(locations.blastDir, f"concatenated{sequenceSuffix}")
    concatenate_sequences(args.sequenceFiles, concatFileName)
    
    # Establish the MMseqs2 query and target databases
    queryDB = ZS_BlastIO.MM_DB(concatFileName, os.path.dirname(args.mms2Exe), locations.tmpDir,
                               args.sequenceType, args.threads)
    queryDB.generate()
    queryDB.index()
    
    targetDB = ZS_BlastIO.MM_DB(targetLink, os.path.dirname(args.mms2Exe), tmpDir,
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

def smain(args, locations):
    # Set up salmon directory
    os.makedirs(locations.salmonDir, exist_ok=True)
    
    # Locate read files
    forwardReads, reverseReads, sampleNames = locate_read_files(args.readsDir, args.readsSuffix,
                                                               args.singleEnd)
    
    # Concatenate all FASTA files into a single file
    concatFileName = os.path.join(locations.salmonDir, f"concatenated.cds") # map to CDS
    concatenate_sequences(args.sequenceFiles, concatFileName)
    
    # Establish the salmon target database
    targetDB = ZS_MapIO.Salmon_DB(concatFileName, os.path.dirname(args.salmonExe),
                                  args.threads)
    salmonDBName = concatFileName + ".salmonDB"
    if not os.path.exists(salmonDBName) or not os.path.exists(salmonDBName + ".ok"):
        print(f"# Indexing '{concatFileName}' for salmon read quantification...")
        targetDB.index()
        touch_ok(salmonDBName)
    
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
        sampleDir = os.path.join(locations.salmonDir, sampleName)
        os.makedirs(sampleDir, exist_ok=True)
        
        # Run salmon (if not already run)
        if not os.path.exists(sampleDir + ".ok"):
            print(f"# Running salmon quantification for '{sampleName}'...")
            salmon.quant(sampleDir)
            touch_ok(sampleDir + ".ok")
    
    print("Salmon quant complete!")

def fmain(args, locations):
    # Set up filter directory for this run
    os.makedirs(locations.filterDir, exist_ok=True)
    
    filterRunName = os.path.basename(args.runDir)
    filterRunDir = os.path.join(locations.filterDir, filterRunName)
    os.makedirs(filterRunDir, exist_ok=True)
    
    mostRecentDir = os.path.join(locations.filterDir, "most_recent")
    if os.path.exists(mostRecentDir) or os.path.islink(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(filterRunDir, mostRecentDir)
    
    # Figure out if we've already filtered this run and exit if so
    outputFileName = os.path.join(filterRunDir, locations.filteredClusterFile)
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"A filtered BINge output file '{locations.filteredClusterFile}' " +
                              f"and its .ok file already exists within '{filterRunDir}'; " +
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
        touch_ok(outputFileName)
        return # exit the program here
    
    # Load transcripts into memory for quick access
    transcriptRecords = FastaCollection(args.sequenceFiles)
    
    # Parse BLAST results (if relevant)
    if args.useBLAST:
        blastResults = ZS_BlastIO.BLAST_Results(args.blastFile)
        blastResults.evalue = args.evalue
        blastResults.num_hits = 1 # only need to keep the best hit for each sequence
        blastResults.parse()
        blastDict = blastResults.results
    
    # Parse annotation file (if relevant)
    if args.useGFF3:
        annotIDs = set()
        for gff3File in args.gff3Files:
            annotIDs = annotIDs.union(parse_gff3_ids(gff3File))
    else:
        annotIDs = set()
    
    # Parse salmon quant (if relevant)
    if args.useSalmon:
        sampleNames = [ f"{i}" for i in range(len(args.salmonFiles))] # sample names don't matter
        quantCollection = parse_quants(args.salmonFiles, sampleNames)
    
        # Get an appropriate read count value to use as a cut-off
        """Note: we want this to be somewhat strict since each cut-off SAVES a sequence,
        so if this is the only thing that passes on a cluster, we want to be confident
        that it's actually good."""
        countCutoff = get_counts_cutoff_by_percentiles(
            binnedDict, unbinnedDict,
            quantCollection
        )
    
    # Perform filtration of clusters with available evidence
    toDropBinned = set()
    toDropUnbinned = set()
    for index, toFilterDict in enumerate([binnedDict, unbinnedDict]):
        for clusterID, seqIDs in toFilterDict.items():
            # Skip unless we are filtering binned clusters
            if index == 0 and (not args.filterBinned):
                continue
            
            # Check 1: Retain if it has a reference sequence
            if args.useGFF3:
                if any([ seqID in annotIDs for seqID in seqIDs ]):
                    continue
            
            # Handle read alignment values
            if args.useSalmon:
                transcriptCounts = [
                    quantCollection.get_transcript_count(seqID)
                    for seqID in seqIDs
                    if seqID in quantCollection.quant
                ]
                
                # Check 2: FILTER if it lacks 1x coverage
                has1x = determine_if_1x_filter(transcriptCounts, seqIDs,
                                               transcriptRecords, args.readLength)
                
                if not has1x:
                    if index == 0:
                        toDropBinned.add(clusterID)
                    else:
                        toDropUnbinned.add(clusterID)
                    continue # since we're filtering, we toDrop it then continue
                
                # Check 3: Retain if any have good read alignment in >1 sample
                if len(quantCollection.samples) > 1:
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
            if args.useBLAST:
                if any([ seqID in blastDict for seqID in seqIDs ]): # dict is evalue filtered already
                    continue
            
            # Check 5: Retain if the sequence is long enough
            "Note: we include stop codons in the length calculations here"
            if any(
            [
                ( len(str(transcriptRecords[seqID])) / 3 ) >= args.minimumLength
                for seqID in seqIDs
            ]):
                continue
            
            # If we make it here, this cluster should be filtered
            if index == 0:
                toDropBinned.add(clusterID)
            else:
                toDropUnbinned.add(clusterID)
    
    # Merge cluster dicts back together
    toDrop = toDropBinned.union(toDropUnbinned)
    clusterDict = merge_cluster_dicts(binnedDict, unbinnedDict, toDrop)
    
    # Write filtered clusters to file
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines
        for clusterNum, value in clusterDict.items():
            clusterIDs, binType = value
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum}\t{seqID}\t{binType}\n")
        touch_ok(outputFileName)
    
    # Print some statistics for the user
    print("# BINge_post.py 'filter' statistics:")
    if args.useSalmon:
        print(f"> A salmon quant value >= {countCutoff} resulted in cluster retention")
    print(f"> Removed {len(toDropBinned)} binned clusters")
    print(f"> Removed {len(toDropUnbinned)} unbinned clusters")
    print(f"> Result contains {len(clusterDict)} clusters")
    
    print("Filtering complete!")

def rmain(args, locations):
    # Set up representatives directory for this run
    os.makedirs(locations.representativesDir, exist_ok=True)
    
    reprRunName = os.path.basename(os.path.dirname(args.bingeFile))
    reprRunDir = os.path.join(locations.representativesDir, reprRunName)
    os.makedirs(reprRunDir, exist_ok=True)
    
    mostRecentDir = os.path.join(locations.representativesDir, "most_recent")
    if os.path.exists(mostRecentDir) or os.path.islink(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(reprRunDir, mostRecentDir)
    
    # Figure out if we've already gotten representatives for this run and exit if so
    for fileName in [locations.representativeMRNA, locations.representativeCDS, locations.representativeAA]:
        outputFileName = os.path.join(reprRunDir, fileName)
        if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
            raise FileExistsError(f"A BINge representative file '{fileName}' " +
                                  f"and its .ok file already exists within '{reprRunDir}'; " +
                                  "this program will not overwrite an existing file. " +
                                  "To resume or overwrite, move/rename/delete this file then try again.")
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(args.bingeFile, "all")
    
    # Load transcripts into memory for quick access
    mrnaRecords = FastaCollection(args.mrnaSequenceFiles)
    cdsRecords = FastaCollection(args.cdsSequenceFiles)
    aaRecords = FastaCollection(args.aaSequenceFiles)
    
    # Parse BLAST results (if relevant)
    if args.useBLAST:
        blastResults = ZS_BlastIO.BLAST_Results(args.blastFile)
        blastResults.evalue = args.evalue
        blastResults.num_hits = 1 # only need to keep the best hit for each sequence
        blastResults.parse()
        blastDict = blastResults.results
    else:
        blastDict = {}
    
    # Parse annotation file (if relevant)
    if args.useGFF3:
        annotIDs = set()
        for gff3File in args.gff3Files:
            annotIDs = annotIDs.union(parse_gff3_ids(gff3File))
    else:
        annotIDs = set()
    
    # Parse salmon quant (if relevant)
    if args.salmonFiles != []:
        sampleNames = [ f"{i}" for i in range(len(args.salmonFiles))] # sample names don't matter
        quantCollection = parse_quants(args.salmonFiles, sampleNames)
    else:
        quantCollection = None
    
    # Write cluster representatives to file
    with open(os.path.join(reprRunDir, locations.representativeMRNA), "w") as mrnaOut, \
    open(os.path.join(reprRunDir, locations.representativeCDS), "w") as cdsOut, \
    open(os.path.join(reprRunDir, locations.representativeAA), "w") as aaOut:
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
                    # Raise error for non-existing sequences
                    if not seqID in cdsRecords:
                        raise KeyError(f"'{seqID}' not found in any input FASTA files; have you modified " +
                                       "any files after initialisation or clustering?")
                    
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
                            blastDict[seqID][0][7] if seqID in blastDict else 0, # best hit is index 0
                            counts,
                            len(str(cdsRecords[seqID])),
                            seqID
                        ]
                    else:
                        thisEvidenceList = [
                            0,
                            blastDict[seqID][0][7] if seqID in blastDict else 0, # bitscore is index 7
                            counts,
                            len(str(cdsRecords[seqID])),
                            seqID
                        ]
                    evidenceLists.append(thisEvidenceList)
                
                # Sort our evidence lists in a way where the first value is the best
                evidenceLists.sort(
                    key = lambda x: (
                        -x[0], -x[1], -x[2], -x[3]
                    )
                )
                
                # Choose the top sequence
                seqID = evidenceLists[0][-1]
            
            # Write to file
            ## mRNA
            mrnaSeq = str(mrnaRecords[seqID])
            mrnaString = format_representative(clusterID, seqID, mrnaSeq)
            mrnaOut.write(mrnaString)
            
            ## CDS
            cdsSeq = str(cdsRecords[seqID])
            cdsString = format_representative(clusterID, seqID, cdsSeq)
            cdsOut.write(cdsString)
            
            ## AA
            aaSeq = str(aaRecords[seqID])
            aaString = format_representative(clusterID, seqID, aaSeq)
            aaOut.write(aaString)
        
        # Touch .ok files after writing all outputs
        touch_ok(os.path.join(reprRunDir, locations.representativeMRNA))
        touch_ok(os.path.join(reprRunDir, locations.representativeCDS))
        touch_ok(os.path.join(reprRunDir, locations.representativeAA))
    
    print("Representative picking complete!")

def dmain(args, locations):
    # Set up DGE directory for this run
    os.makedirs(locations.dgeDir, exist_ok=True)
    
    dgeRunName = os.path.basename(os.path.dirname(args.bingeFile))
    dgeRunName = os.path.join(locations.dgeDir, dgeRunName)
    os.makedirs(dgeRunName, exist_ok=True)
    
    mostRecentDir = os.path.join(locations.dgeDir, "most_recent")
    if os.path.exists(mostRecentDir) or os.path.islink(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(dgeRunName, mostRecentDir)
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(args.bingeFile, "all")
    
    # Generate tx2gene file (if not already done)
    tx2geneFile = os.path.join(dgeRunName, locations.tx2geneFile)
    if not os.path.exists(tx2geneFile) or not os.path.exists(tx2geneFile + ".ok"):
        print(f"# Generating '{locations.tx2geneFile}' file for cluster count summarisation...")
        write_tx2gene(tx2geneFile, clusterDict)
        touch_ok(tx2geneFile)
    
    # Derive the directories from the salmon files
    salmonDirs = [ os.path.dirname(salmonFile) for salmonFile in args.salmonFiles ]
    
    # Tabulate salmon QC metrics (if not already done)
    salmonQCFile = os.path.join(dgeRunName, locations.salmonQCFile)
    if not os.path.exists(salmonQCFile) or not os.path.exists(salmonQCFile + ".ok"):
        print(f"# Generating '{locations.salmonQCFile}' file for salmon QC assessment...")
        write_salmonQC(salmonDirs, salmonQCFile, clusterDict)
        touch_ok(salmonQCFile)
    
    # Write list of samples for DGE analysis (if not already done)
    sampleFile = os.path.join(dgeRunName, locations.salmonSampleFile)
    if not os.path.exists(sampleFile) or not os.path.exists(sampleFile + ".ok"):
        print(f"# Generating '{locations.salmonSampleFile}' file for sample name listing...")
        write_samples(salmonDirs, sampleFile)
        touch_ok(sampleFile)
    
    # Generate R script for DGE analysis (if not already done)
    rScriptFile = os.path.join(dgeRunName, locations.rScriptFile)
    if not os.path.exists(rScriptFile) or not os.path.exists(rScriptFile + ".ok"):
        print(f"# Generating '{locations.rScriptFile}' file for DGE analysis...")
        write_r_script(locations.salmonDir, sampleFile, tx2geneFile, rScriptFile)
        touch_ok(rScriptFile)
    
    print("DGE preparation complete!")

def amain(args, locations):
    # Set up annotate directory for this run
    os.makedirs(locations.annotateDir, exist_ok=True)
    
    annotateRunName = os.path.basename(os.path.dirname(args.bingeFile))
    annotateRunDir = os.path.join(locations.dgeDir, annotateRunDir)
    os.makedirs(annotateRunDir, exist_ok=True)
    
    mostRecentDir = os.path.join(locations.annotateDir, "most_recent")
    if os.path.exists(mostRecentDir) or os.path.islink(mostRecentDir):
        os.unlink(mostRecentDir)
    os.symlink(annotateRunDir, mostRecentDir)
    
    # Figure out if we've already annotated this run and exit if so
    outputFileName = os.path.join(annotateRunDir, locations.annotationFile)
    if os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok"):
        raise FileExistsError(f"A BINge annotation file '{locations.annotationFile}' " +
                              f"and its .ok file already exists within '{annotateRunDir}'; " +
                              "this program will not overwrite an existing file. " +
                              "To resume or overwrite, move/rename/delete this file then try again.")
    
    # Validate representatives file
    "Just use the .aa file since it will be the smallest and IDs are the same across all files"
    representativesFasta = os.path.join(locations.representativesDir, annotateRunName,
                                        locations.representativeAA)
    if not os.path.exists(representativesFasta) or not os.path.exists(representativesFasta + ".ok"):
        raise FileNotFoundError(f"Unable to locate '{representativesFasta}' or '{representativesFasta}.ok'; " +
                                "have you run 'representatives' yet?")
    
    # Parse the representatives file to find the sequences we need
    clustToRep, repToClust = parse_binge_representatives(representativesFasta)
    
    # Parse .obo file
    "Parse early to error out if format is incorrect"
    goObo = obo_parser.GODag(args.oboFile)
    
    # Step 1: initial parse of BLAST file to begin formatting the annotation table
    step1File = os.path.join(annotateRunDir, "tmp.step1.tsv") # will overwrite if exists
    hitMapDict = init_table(
        clustToRep, repToClust,
        args.blastFile, args.evalue, args.numHits,
        step1File, args.databaseTag, args.largeTable)
    
    # Step 2: parse idmapping file
    parse_idmap(args.idmappingFile, hitMapDict)
    
    # Step 3: update the annotation table with GOs
    step3File = os.path.join(annotateRunDir,  "tmp.step3.tsv")
    update_table_with_gos(step1File, step3File, hitMapDict, goObo)
    os.unlink(step1File) # clean up first temporary file now that it has been used
    
    # Step 4: update the annotation table a final time to include sequence details
    update_table_with_seq_details(step3File, outputFileName, hitMapDict,
                                  representativesFasta, args.targetFile)
    os.unlink(step3File) # now clean up the second temporary file
    touch_ok(outputFileName)
    
    print("Annotation complete!")

if __name__ == "__main__":
    main()
