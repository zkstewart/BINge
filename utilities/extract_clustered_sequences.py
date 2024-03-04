#! python3
# extract_clustered_sequences.py

# Script to take in the BINge working directory (which should have .nucl 
# files corresponding to mRNA/CDS sequences in it) and generate a single
# FASTA file containing the sequences that were clustered.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.fasta_handling import FastaCollection
from modules.parsing import parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input locations
    if not os.path.isdir(args.bingeDir):
        print(f'I am unable to locate the BINge working directory ({args.bingeDir})')
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
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
    usage = """%(prog)s will take the BINge clustering result file, which may have been
    filtered, alongside the directory where BINge was initially run. It will extract ALL
    sequences present in your cluster file and write them to a single FASTA file. This
    is intended to be useful for performing read mapping against the clustered sequences.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input BINge clustering result file")
    p.add_argument("-d", dest="bingeDir",
                   required=True,
                   help="Input location where BINge was run")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for sequence FASTA")
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeResultFile,
        "all"
    )
    
    # Parse the .nucl files in the BINge working directory
    transcriptRecords = FastaCollection([
        os.path.join(args.bingeDir, f)
        for f in os.listdir(args.bingeDir)
        if f.endswith(".nucl")
    ])
    
    # Extract the clustered sequences
    try:
        with open(args.outputFileName, "w") as fileOut:
            for clusterIDs in clusterDict.values():
                for seqID in clusterIDs:
                    try:
                        record = transcriptRecords[seqID]
                    except KeyError as e:
                        print(f"ERROR: Unable to locate sequence '{seqID}' in the BINge working directory.")
                        raise e # Raise the error again for clean up
                    
                    fileOut.write(f">{seqID}\n{str(record)}\n")
    except Exception as e:
        print("ERROR: Something went wrong, so the output file is being removed.")
        
        # Clean up the output file
        if os.path.isfile(args.outputFileName):
            os.unlink(args.outputFileName)
        
        # Provide the error message if applicable
        if type(e) == KeyError:
            pass # we don't need to raise the error since we got a useful error msg already
        else:
            print(f"ERROR: Exception message was: {e}")
        quit()
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
