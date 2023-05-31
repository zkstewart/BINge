#! python3
# BINge.py
# BIN Genes for Expression analyses - filter module

# Utility program to follow downstream of BINge clustering, allowing one to
# filter clusters to limit results only to sequences with good support.

import os, argparse, sys
import numpy as np
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from BINge_counter import validate_salmon_files, parse_binge_clusters, \
    parse_equivalence_classes, parse_quants
from BINge_tuning import validate_cluster_file
from BINge_representatives import FastaCollection

from Various_scripts import ZS_BlastIO, ZS_ClustIO, ZS_SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeFile):
        print(f'I am unable to locate the BINge cluster file ({args.bingeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for fastaFile in args.fastaFiles:
        if not os.path.isfile(fastaFile):
            print(f'I am unable to locate the transcript FASTA file ({fastaFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate cluster file
    args.isBinge = validate_cluster_file(args.bingeFile)
    
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
    
    # Validate evalue and its logical necessity
    if args.blastFile == None and "--evalue" in sys.argv:
        print("You've used --evalue despite not giving a BLAST file as input.")
        print("Just in case you've forgotten to provide the BLAST file, I'm going to exit now.")
        print("Either drop the --evalue parameter, or include a --blast parameter.")
        quit()
    if args.evalue < 0:
        print("--evalue must be a positive float value (minimum value of zero)")
        print("Make sure to fix this and try again.")
        quit()
    
    # Validate salmon files if relevant
    if args.salmonFiles != []:
        args.salmonFileFormat = validate_salmon_files(args.salmonFiles)
    
    # Validate length cut-off
    if args.minimumLength < 0:
        print("--length must be an integer with a minimum value of zero")
        print("Make sure to fix this and try again.")
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
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

def get_bottom_percentile_of_counts(clusterDictList, salmonCollection, salmonFileFormat,
                                    percentile=50, minCount=0.5):
    '''
    Parses the sequence IDs out of one or more cluster dictionaries, alongside a salmon
    EquivalenceClassCollection or QuantCollection (with format indicated), and calculates
    the nth percentile of read count values in the dataset.
    
    Parameters:
        clusterDictList -- a list of one or more cluster dictionaries with format like:
                           {
                               1: [ seqid1, seqid2, ... ],
                               2: [...],
                               ...
                           }
        salmonCollection -- an EquivalenceClassCollection or QuantCollection.
        salmonFileFormat -- a string == "ec" if the collection is an EquivalenceClassCollection,
                            and "quant" if a QuantCollection
        percentile -- an integer for use with np.percentile.
        minCount -- a float specifying the minimum value we'll use if the percentile gives
                    a read count of zero (0).
    '''
    counts = [
        np.mean(salmonCollection.get_transcript_count(seqID))
        
        for clusterDict in clusterDictList
        for seqIDs in clusterDict.values()
        for seqID in seqIDs
        if (salmonFileFormat == "ec" and seqID in salmonCollection.ids)
        or (salmonFileFormat == "quant" and seqID in salmonCollection.quant)
    ]
    fiftyPercentile = np.percentile(counts, percentile)
    countCutoff = fiftyPercentile if fiftyPercentile > 0 else minCount # need at least SOME alignment
    
    return countCutoff

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

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside several other files to filter predicted
    clusters, retaining only those which meet specific criteria defined herein.
    
    Optional inclusions for cluster filtering include:
    1) A BLAST outfmt6 file to include best-BLAST evidence for cluster filtering.
    2) A GFF3 of the reference organism used as part of BINge clustering, or just a text
    file listing the IDs of its sequences.
    3) One or more Salmon quant.sf or equivalence class files indicating the read alignment
    made to transcripts.
    
    Clusters will be retained if they meet any of the below criteria in this order of
    importance i.e., if they contain sequence(s):
    1) from the reference GFF3, or
    2) with a significant BLAST E-value, or
    3) with read alignment > cut-off value from more than 1 sample (if multiple
    Salmon files are given), or
    4) have an ORF >= cut-off length.
    
    As an example, if filter 2 passes i.e., a sequence in a cluster has a significant
    BLAST hit, that cluster will NOT be filtered even if its ORF length is less than
    that cut-off.
    
    Hence, without providing any optional files, it just filters clusters to those with
    a minimum sequence length. Ideally, you should provide all the optional files to
    make filtering decisions based on biological evidence. 
    
    The output is modified BINge cluster TSV file with failing clusters removed. The cluster
    numbers will be revised to account for absences, so please expect the cluster numbers to
    change after filtering.
    
    Note: for testing purposes, this script will accept either a BINge or CD-HIT cluster
    file as input (-i). Also note that the input file can be in any format (nucleotide mRNA,
    CDS, or protein translation) but this code will work best if you give the CDS or protein
    sequence since that makes it easier to get the correct ORF length. Otherwise, ORFs will
    be predicted in a very simplistic way.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeFile",
                   required=True,
                   help="Input BINge cluster file")
    p.add_argument("-f", dest="fastaFiles",
                   required=True,
                   nargs="+",
                   help="Input one or more FASTAs containing transcripts listed in the cluster file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    p.add_argument("-s", dest="sequenceFormat",
                   required=True,
                   choices=["transcript", "cds", "protein"],
                   help="""Indicate whether your file contains transcripts (i.e., inclusive
                   of UTR; as nucleotide), CDS (as nucleotide), or protein sequences.""")
    # Optional
    p.add_argument("--blast", dest="blastFile",
                   required=False,
                   help="An outfmt6 file from BLAST or MMseqs2 against a relevant database.",
                   default=None)
    p.add_argument("--annot", dest="annotationFile",
                   required=False,
                   help="A GFF3 or text file containing the representative organism's sequence IDs",
                   default=None)
    p.add_argument("--salmon", dest="salmonFiles",
                   nargs="+",
                   required=False,
                   help="One or more salmon eq_classes.txt or quant.sf files (one or the other)",
                   default=[])
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="""If you have provided a BLAST file, optionally specify an E-value
                   threshold to use when considering a hit to be statistically significant.""",
                   default=1e-5)
    p.add_argument("--length", dest="minimumLength",
                   required=False,
                   type=int,
                   help="""Specify the minimum length of an ORF (in amino acids) that must
                   be met for a cluster to pass filtration (if it hasn't passed any of the
                   previous filter checks already).""",
                   default=50)
    p.add_argument("--filter_binned", dest="filterBinned",
                   required=False,
                   action="store_true",
                   help="""Optionally, extend filtering to clusters that were binned. By default,
                   binned clusters are assumed to be good since they've already met sequence
                   similarity criteria to the reference genome.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the CD-HIT / BINge cluster file
    if args.isBinge:
        binnedDict = parse_binge_clusters(args.bingeFile, "binned")
        unbinnedDict = parse_binge_clusters(args.bingeFile, "unbinned")
    else:
        binnedDict = {} # CD-HIT doesn't give us binned/unbinned information
        unbinnedDict = ZS_ClustIO.CDHIT.parse_clstr_file(args.bingeFile)
    
    # Load transcripts into memory for quick access
    transcriptRecords = FastaCollection(args.fastaFiles)
    
    # Parse BLAST results (if relevant)
    if args.blastFile != None:
        blastResults = ZS_BlastIO.BLAST_Results(args.blastFile)
        blastResults.evalue = args.evalue
        blastResults.num_hits = 1 # only need to keep the best hit for each sequence
        blastResults.parse_blast_hit_coords()
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
    
    # Parse salmon counts (if relevant)
    if args.salmonFiles != []:
        sampleNames = [ f"{i}" for i in range(len(args.salmonFiles))] # sample names don't matter
        
        if args.salmonFileFormat == "ec":
            salmonCollection = parse_equivalence_classes(args.salmonFiles, sampleNames)
        else:
            salmonCollection = parse_quants(args.salmonFiles, sampleNames)
    else:
        salmonCollection = None
    
    # Get the bottom percentile of read counts to use as a cut-off
    """Note: we want this to be somewhat strict since each cut-off SAVES a sequence,
    so if this is the only thing that passes on a cluster, we want to be confident
    that it's actually good. The 50th percentile seems pretty good for deciding this."""
    countCutoff = get_bottom_percentile_of_counts(
        [binnedDict, unbinnedDict],
        salmonCollection,
        args.salmonFileFormat,
        50, 0.5
    )
    
    # Perform filtration of clusters with available evidence
    toDrop = set()
    if not args.filterBinned:
        toFilterDicts = [unbinnedDict]
    else:
        toFilterDicts = [binnedDict, unbinnedDict]
    
    for toFilterDict in toFilterDicts:
        for clusterID, seqIDs in toFilterDict.items():
            # Check 1: Retain if it has a reference sequence
            if any([ seqID in annotIDs for seqID in seqIDs ]):
                continue
            
            # Check 2: Retain if any have a significant E-value
            if any([ seqID in blastDict for seqID in seqIDs ]): # if in blastDict, it has a good E-value
                continue
            
            # Check 3: Retain if any have good read alignment in >1 sample
            if len(salmonCollection.samples) > 1:
                if any(
                [
                    sum([ countValue > countCutoff for countValue in salmonCollection.get_transcript_count(seqID) ]) > 1
                    for seqID in seqIDs
                    if (args.salmonFileFormat == "ec" and seqID in salmonCollection.ids)
                    or (args.salmonFileFormat == "quant" and seqID in salmonCollection.quant)
                ]):
                    continue
            
            # Check 3 (alt): Retain if it has good read alignment in the 1 sample
            else:
                if any(
                [
                    salmonCollection.get_transcript_count(seqID)[0] > countCutoff
                    for seqID in seqIDs
                    if (args.salmonFileFormat == "ec" and seqID in salmonCollection.ids)
                    or (args.salmonFileFormat == "quant" and seqID in salmonCollection.quant)
                ]):
                    continue
            
            # Check 4: Retain if the sequence is long enough
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
    print(f"BINge_filter removed {len(toDrop)} clusters")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
