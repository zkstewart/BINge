#! python3
# BINge_filter.py
# BIN Genes for Expression analyses - filter module

# Utility program to follow downstream of BINge clustering, allowing one to
# filter clusters to limit results only to sequences with good support.

import os, argparse, sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts.Function_packages import ZS_BlastIO

from modules.fasta_handling import ZS_SeqIO, FastaCollection
from modules.validation import validate_salmon_files, validate_cluster_file
from modules.parsing import parse_equivalence_classes, parse_quants, parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeFile):
        print(f'I am unable to locate the BINge cluster file ({args.bingeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for filterFile in args.filterFiles:
        if not os.path.isfile(filterFile) and not os.path.islink(filterFile):
            print(f'I am unable to locate the filter file ({filterFile)')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
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
    
    # Validate cluster file
    args.isBinge = validate_cluster_file(args.bingeFile)
    if not args.isBinge:
        print("CD-HIT is no longer supported by BINge_filter.py")
        print("Sorry about that - make sure to use BINge results herein.")
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
    
    # Validate readLength and its logical necessity
    if args.require1x == False and "--readLength" in sys.argv:
        print("You've used --readLength despite not giving the --require1x flag.")
        print("Just in case you're confused about whether you want this filtration to occur " +
              "or not, I'm going to exit now.")
        print("Either drop the --readLength parameter, or include --require1x.")
        quit()
    if args.readLength < 1:
        print("--readLength must be a positive integer value (minimum value of 1)")
        print("Make sure to fix this and try again.")
        quit()
    if args.readLength < 50:
        print(f"--readLength has received an unusually low value ({args.readLength})")
        print("We will keep running, but this may result in excessive filtration.")
    
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
                    fastaIDs.add(line.rstrip("\r\n ").split(" ")[0])
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

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside several other files to filter predicted
    clusters, retaining only those which meet specific criteria defined herein.
    
    You should provide the -f argument the location of the BINge results directory since
    it will contain all of the annotation.nucl and transcriptome.nucl files necessary for
    this script to work.
    
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
    
    Note: the input file can be in any format (nucleotide mRNA, CDS, or protein translation)
    but this code will work best if you give the CDS or protein sequence since that makes it
    easier to get the correct ORF length. Otherwise, ORFs will be predicted in a very
    simplistic way.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeFile",
                   required=True,
                   help="Input BINge cluster file")
    p.add_argument("-f", dest="fastaDir",
                   required=True,
                   help="Input directory containing FASTAs listed in the cluster file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    p.add_argument("-s", dest="sequenceFormat",
                   required=True,
                   choices=["transcript", "cds", "protein"],
                   help="""Indicate whether your file contains transcripts (i.e., inclusive
                   of UTR; as nucleotide), CDS (as nucleotide), or protein sequences.""")
    # Optional
    p.add_argument("--fastas_filter", dest="filterFiles",
                   nargs="+",
                   required=False,
                   help="""Optionally, specify one or more FASTA file(s) to retain clusters
                   ONLY if they contain any of the sequences in these file(s); this filter 
                   supercedes all others.""",
                   default=[])
    p.add_argument("--annot", dest="annotationFile",
                   required=False,
                   help="""Optionally, specify a GFF3 or text file containing the representative
                   organism's sequence IDs to retain an unbinned cluster if it contains any of
                   these sequences.""",
                   default=None)
    p.add_argument("--length", dest="minimumLength",
                   required=False,
                   type=int,
                   help="""Specify the minimum length of an ORF (in amino acids) that must
                   be met for a cluster to pass filtration (if it hasn't passed any other
                   filter checks already); make this strict! (default==150)""",
                   default=150)
    p.add_argument("--blast", dest="blastFile",
                   required=False,
                   help="""Optionally, specify an outfmt6 file from BLAST or MMseqs2 against a
                   relevant database to use for filtering purposes.""",
                   default=None)
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="""If you have provided a BLAST file, optionally specify an E-value
                   threshold to use when considering a hit to be statistically significant
                   (default==1e-10).""",
                   default=1e-10)
    p.add_argument("--salmon", dest="salmonFiles",
                   nargs="+",
                   required=False,
                   help="""Optionally, specify one or more salmon eq_classes.txt or quant.sf
                   files (one or the other type only) to use with read depth retention""",
                   default=[])
    p.add_argument("--require1x", dest="require1x",
                   required=False,
                   action="store_true",
                   help="""Optionally, require a transcript in a cluster to have at least 1x
                   coverage summed across samples to retain the cluster.""",
                   default=False)
    p.add_argument("--readLength", dest="readLength",
                   required=False,
                   type=int,
                   help="""If specifying --require1x, indicate the length of your reads in
                   basepairs (default==150; just consider the single end length)""",
                   default=150)
    p.add_argument("--filter_binned", dest="filterBinned",
                   required=False,
                   action="store_true",
                   help="""Optionally, extend filtering to clusters that were binned. By default,
                   binned clusters are assumed to be good since they've already met sequence
                   similarity criteria to the reference genome.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    binnedDict = parse_binge_clusters(args.bingeFile, "binned")
    unbinnedDict = parse_binge_clusters(args.bingeFile, "unbinned")
    
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
    if not args.filterBinned:
        toFilterDicts = [unbinnedDict]
    else:
        toFilterDicts = [binnedDict, unbinnedDict]
    
    for toFilterDict in toFilterDicts:
        for clusterID, seqIDs in toFilterDict.items():
            # Check 0: FILTER if it does not contain a filter ID
            if args.filterFiles != []:
                if not any([ seqID in filterIDs for seqID in seqIDs ]):
                    toDrop.add(clusterID)
                    continue # since we're filtering, we toDrop it then continue
            
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
    print("# BINge_filter: Filtration statistics")
    print(f" > used 1x criteria to filter clusters = {args.require1x}")
    print(f" > number of reads = {countCutoff} required for read count rescue")
    print(f" > removed {len(toDrop)} clusters")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
