#! python3
# tx2gene.py

# Utility program for converting a BINge cluster file into
# the tx2gene mapping format used by DESeq2 when summing
# transcript abundance estimates from Salmon/others to the gene level.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.parsing import parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeResultFile):
        print(f'I am unable to locate the BINge result file ({args.bingeResultFile})')
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
    usage = """%(prog)s provides a simple conversion service, taking a BINge clusters
    TSV and converting it into a tx2gene file that can be used with DESeq2's tximport
    in R to load transcript abundances e.g., from Salmon, and tally them to the gene level.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name")
    # Optional
    p.add_argument("--only_binned", dest="onlyBinned",
                   required=False,
                   action="store_true",
                   help="Optionally, only output tx2gene mapping for binned sequence clusters",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeResultFile,
        "all" if args.onlyBinned == False else "binned"
    )
    
    # Create the tx2gene file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("TXNAME\tGENEID\n")
        for clusterNum, seqIDs in clusterDict.items():
            for seqID in seqIDs:
                fileOut.write(f"{seqID}\tcluster-{clusterNum}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
