#! python3
# BINge_representatives.py
# BIN Genes for Expression analyses - representatives module

# Utility program to select a representative sequence for each cluster
# predicted by BINge making use of (potentially) several lines of evidence.

import os, argparse

from modules.fasta_handling import FastaCollection
from modules.validation import validate_salmon_files, validate_cluster_file
from modules.parsing import parse_equivalence_classes, parse_quants, parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeFile):
        print(f'I am unable to locate the BINge cluster file ({args.bingeFile})')
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
        print("CD-HIT is no longer supported by BINge_representatives.py")
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
    
    # Validate salmon files if relevant
    if args.salmonFiles != []:
        args.salmonFileFormat = validate_salmon_files(args.salmonFiles)
    
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

def parse_blast_outfmt6(blastFile, evalueCutoff):
    '''
    Simply parses a BLAST outfmt6 file to get the best E-value hit for each query sequence.
    
    Parameters:
        blastFile -- a string indicating the location of a BLAST results file in
                     outfmt6 format.
    Returns:
        blastDict -- a dict with structure like:
                     {
                         'seqID1': evalue,
                         'seqID2': evalue,
                         ...
                     }
    '''
    blastDict = {}
    
    with open(blastFile, "r") as fileIn:
        for line in fileIn:
            # Extract details
            sl = line.rstrip("\r\n").split('\t')
            queryID = sl[0]
            evalue = float(sl[10])
            
            # Skip if evalue isn't significant
            if evalue > evalueCutoff:
                continue
            
            # Store result if none exist
            if queryID not in blastDict:
                blastDict[queryID] = evalue
            # Overwrite result if this is better
            elif evalue < blastDict[queryID]:
                blastDict[queryID] = evalue
    
    return blastDict

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

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside several other files to pick out a single
    representative sequence from each cluster. These representative sequences can be used
    for downstream functions including gene annotation.
    
    You should provide the -f argument the location of the BINge results directory since
    it will contain all of the annotation.nucl and transcriptome.nucl files necessary for
    this script to work.
    
    Optional inclusions for representative picking include:
    1) A BLAST outfmt6 file to select a representative including best-BLAST evidence.
    2) A GFF3 of the reference organism used as part of BINge clustering, or just a text
    file listing the IDs of its sequences.
    3) One or more Salmon quant.sf or equivalence class files indicating the read alignment
    made to transcripts.
    
    Evidence weighting for picking a representative is as such:
    (is a reference organism sequence) > (BLAST E-value) > (number of reads aligned) >
    sequence length.
    
    Hence, without providing any optional files, it just picks a
    representative on the basis of sequence length (longest being best). Ideally, you
    should provide at least the BLAST and/or Salmon files to pick a representative using
    more relevant biological criteria.
    
    The output is a FASTA file containing these representatives. The ID for each sequence
    follows a format like '>Cluster-1 representative=transcript_92'.
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
                   default=1e-3)
    p.add_argument("--only_binned", dest="onlyBinned",
                   required=False,
                   action="store_true",
                   help="Optionally, only output representatives for binned sequence clusters",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeFile,
        "all" if args.onlyBinned == False else "binned"
    )
    
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
                            blastDict[seqID] if seqID in blastDict else 1, # stores E-value
                            counts,
                            len(str(transcriptRecords[seqID])),
                            seqID
                        ]
                    else:
                        thisEvidenceList = [
                            0,
                            blastDict[seqID] if seqID in blastDict else 1, # stores E-value
                            counts,
                            len(str(transcriptRecords[seqID])),
                            seqID
                        ]
                    evidenceLists.append(thisEvidenceList)
                
                # Sort our evidence lists in a way where the first value is the best
                evidenceLists.sort(
                    key = lambda x: (
                        -x[0], x[1], -x[2], -x[3] # with E-value, smaller is better so don't '-' it
                    )
                )
                
                # Choose the top sequence
                seqID = evidenceLists[0][-1]
            
            # Write to file
            seq = str(transcriptRecords[seqID])
            fastaString = format_representative(clusterID, seqID, seq)
            fileOut.write(fastaString)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
