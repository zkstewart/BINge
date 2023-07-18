#! python3
# update_outfmt6.py

# Script to take in a BINge_representatives FASTA file alongside the outfmt6
# BLAST results that were used for selecting representatives. It will produce a
# new outfmt6 file where non-representative sequences are removed, and the
# representative sequences have their IDs switched with their cluster ID.

import os, argparse
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.representativesFasta):
        print(f'I am unable to locate the representatives FASTA file ({args.representativesFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.blastFile):
        print(f'I am unable to locate the BLAST outfmt6 file ({args.blastFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

## Main
def main():
    # User input
    usage = """%(prog)s reads in the output FASTA of BINge_representatives.py alongside
    an outfmt6 file from BLAST or MMseqs2. The outfmt6 file will be updated to 1) remove
    any sequences that aren't the representatives of their cluster, and 2) update
    the result IDs to indicate the cluster number (e.g., 'cluster-42') rather than their
    original sequence ID.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-r", dest="representativesFasta",
                   required=True,
                   help="Input representatives FASTA file")
    p.add_argument("-b", dest="blastFile",
                   required=True,
                   help="Input BLAST outfmt6 file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output name for updated outfmt6 file")
    
    args = p.parse_args()
    validate_args(args)
    
    #### Disable functionality of this code until a better solution is found ####
    
    print("update_outfmt6.py has a flaw in its logic")
    print("If a representative is from a reference genome, it may not have a hit")
    print("This would result in NO hit showing through in the results")
    print("You're better off running BLAST/MMSeqs2 again using the output representative sequences" +
          " to prevent any silent errors in your data")
    quit()
    
    #### ----- ####
    
    # Associate representatives with their cluster IDs
    clustersDict = {} # don't like this variable name anymore but I need to keep it ...
    with open(args.representativesFasta, "r") as fileIn: # ... consistent with unify_representatives
        fastaRecords = SeqIO.parse(fileIn, "fasta")
        for record in fastaRecords:
            clusterID = record.id
            representativeID = record.description.split("representative=")[1].split(" ")[0]
            clustersDict[representativeID] = clusterID
    
    # Parse the outfmt6 file and write revised output
    with open(args.blastFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Get data from this line
            sl = line.rstrip("\r\n ").split("\t")
            seqID = sl[0]
            
            # If this line is a representative, format it for output
            if seqID in clustersDict:
                clusterID = clustersDict[seqID]
                fileOut.write("{0}\n".format("\t".join([
                    clusterID,
                    *sl[1:]
                ])))
            # Otherwise, skip over this line
            else:
                continue
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
