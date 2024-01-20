#! python3
# unify_representatives.py

# Script to take in the two classes of files (e.g., CDS and mRNA transcripts)
# that weren't processed with BINge_representatives and unify them to have the 
# same sequences, same IDs, in the same order.

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from BINge_representatives import FastaCollection

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.representativesFasta):
        print(f'I am unable to locate the representatives FASTA file ({args.representativesFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for fileArg in [args.files1, args.files2]:
        for file in fileArg:
            if not os.path.isfile(file):
                print(f"I am unable to locate file at '{file}'")
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
    # Validate output file location
    suffixes = [ f.split(".")[-1].rstrip("._") for f in [args.files1[0], args.files2[0]] ]
    args.outputFileNames = []
    for suffix in suffixes:
        outFileName = args.outputPrefix + "." + suffix
        if os.path.isfile(outFileName):
            print(f'File already exists at output location ({outFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()
        args.outputFileNames.append(outFileName)

## Main
def main():
    # User input
    usage = """%(prog)s reads in the output FASTA of BINge_representatives.py which
    may be in transcript, CDS, or protein format. Here, you can provide one or more
    files from the other two formats that weren't processed with BINge_representatives
    to produce output files with the representative sequences for these other two
    formats.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-r", dest="representativesFasta",
                   required=True,
                   help="Input representatives FASTA file")
    p.add_argument("-f1", dest="files1",
                   required=True,
                   nargs="+",
                   help="Input one or more files of the first format")
    p.add_argument("-f2", dest="files2",
                   required=True,
                   nargs="+",
                   help="Input one or more files of the second format")
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Output prefix for filtered file(s)")
    # Optional
    p.add_argument("--relaxed", dest="beRelaxed",
                   required=False,
                   action="store_true",
                   help="Optionally allow IDs to be missing from one of the input files")
    
    args = p.parse_args()
    validate_args(args)
    
    # Associate representatives with their cluster IDs
    representativesDict = {}
    with open(args.representativesFasta, "r") as fileIn:
        fastaRecords = SeqIO.parse(fileIn, "fasta")
        for record in fastaRecords:
            clusterID = record.id
            representativeID = record.description.split("representative=")[1].split(" ")[0]
            representativesDict[clusterID] = representativeID
    
    # Load transcripts from files1 and files2 into memory for quick access
    files1Records = FastaCollection(args.files1)
    files2Records = FastaCollection(args.files2)
    
    # Write output files
    cleanAndExit = False
    exitedOn = None
    for i in range(2):
        # Exit condition if something went wrong
        if cleanAndExit == True:
            break
        
        # Obtain values for this iteration
        outputFileName = args.outputFileNames[i]
        if i == 0:
            records = files1Records
        else:
            records = files2Records
        
        # Get sequences for each cluster, writing in order
        skippedIDs = []
        with open(outputFileName, "w") as fileOut:
            for clusterID, representativeID in representativesDict.items():
                try:
                    seq = str(records[representativeID])
                except:
                    if not args.beRelaxed:
                        print(f"'{representativeID}' could not be found in your -f{i+1} files.")
                        print("You will need to fix this problem or make sure you're specifying " +
                            "the right files, then try again.")
                        cleanAndExit = True
                        exitedOn = i
                        break
                    else:
                        skippedIDs.append(representativeID)
                        continue
                
                fileOut.write(f">{clusterID} representative={representativeID}\n{seq}\n")
    
    # Clean up if program failed
    if cleanAndExit == True:
        for i in range(0, exitedOn+1):
            outputFileName = args.outputFileNames[i]
            if os.path.isfile(outputFileName):
                os.unlink(outputFileName)
        
        print("Program exited after cleaning up truncated output file(s)!")
    
    # Print IDs if some were skipped
    elif len(skippedIDs) > 0:
        print(f"{len(skippedIDs)} IDs were not found in your input files. These include:")
        for ID in skippedIDs:
            print(ID)
        print("Program otherwise completed successfully!")
    
    # Or just let user know that everything worked well
    else:
        print("Program completed successfully!")

if __name__ == "__main__":
    main()
