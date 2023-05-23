#! python3
# BINge_counter.py
# BIN Genes for Expression analyses

# Utility program for counting reads from Salmon equivalence classes
# and produces a tabulated output.

import os, argparse
from statistics import mean

# Define classes
class EquivalenceClassCollection():
    def __init__(self):
        self.ids = {} # transcript IDs to EC index dict
        self.samples = [] # sample names
        self.ec = {} # equivalence class dicts
        
        self.numTranscripts = None # for checking EC file compatibility
    
    def parse_eq_file(self, eqFile, sample):
        # Validate input parameters
        assert os.path.isfile(eqFile), \
            f"Cannot parse '{eqFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Figure out if we should hold onto transcript IDs or not
        if self.ids == {}:
            needsIDs = True
        else:
            needsIDs = False
        
        # Parse the file proper
        with open(eqFile, "r") as fileIn:
            
            # Parse the first two lines of the file
            numTranscripts = int(fileIn.readline().rstrip("\r\n "))
            numECs = int(fileIn.readline().rstrip("\r\n "))
            
            # Check compatibility of EC files
            if self.numTranscripts == None:
                self.numTranscripts = numTranscripts
            else:
                assert self.numTranscripts == numTranscripts, \
                    "EC files are incompatible since transcript numbers differ!"
            
            # Loop through transcript ID lines
            for i in range(numTranscripts):
                transcriptID = fileIn.readline().rstrip("\r\n ")
                if needsIDs:
                    self.ids[transcriptID] = i
            
            # Loop through EC lines
            thisEC = { i:0 for i in range(0, numTranscripts) } # holds onto counts per transcript for this file
            for _ in range(numECs):
                sl = list(map(int, fileIn.readline().rstrip("\r\n ").split("\t")))
                numIDs, idIndices, numReads = sl[0], sl[1:-1], sl[-1]
                
                # Store relevant data
                for idIndex in idIndices:
                    thisEC[idIndex] += numReads # CAN MAKE TPM COMPATIBLE BY SPLITTING READS HERE
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store EC counts inside parent structure
        for ecIndex, ecCount in thisEC.items():
            self.ec.setdefault(ecIndex, [])
            self.ec[ecIndex].append(ecCount)
    
    def get_transcript_count(self, ecValue, group):
        '''
        Parameters:
            ecValue -- an integer identifying the equivalence class index
                       you want to get counts for, OR a string identifying
                       the transcript ID you want to get counts for
            group -- an integer identifying which group you want to obtain
                     a count for; replicates are automatically averaged and this
                     value is returned
        Returns:
            meanTranscriptCount -- an integer of the mean transcript count
                                   for the indicated replicate group
        '''
        if isinstance(ecValue, int):
            ec = self.ec[ecValue]
        else:
            ec = self.ec[self.ids[ecValue]]
        
        meanTranscriptCount = mean([
            ec[i]
            for i in range(len(self.groups))
            if self.groups[i] == group
        ])
        return meanTranscriptCount
    
    def __repr__(self):
        return "<EquivalenceClassCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeResultFile):
        print(f'I am unable to locate the BINge result file ({args.bingeResultFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for eqFile in args.equivalenceClassFiles:
        if not os.path.isfile(eqFile):
            print(f'I am unable to locate the salmon equivalence class file ({eqFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate the input files are logically sound
    if len(args.equivalenceClassFiles) != len(args.sampleNames):
        print("Your equivalence class and sample names are incompatible!")
        print((f"I'm seeing {len(args.equivalenceClassFiles)} equivalence class files " +
              f"and {len(args.sampleNames)} sample names"))
        print("These numbers should be the same. You need to fix this up and try again.")
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_equivalence_classes(equivalenceClassFiles, sampleNames):
    ecCollection = EquivalenceClassCollection()
    for i in range(len(equivalenceClassFiles)):
        eqFile = equivalenceClassFiles[i]
        sample = sampleNames[i]
        
        ecCollection.parse_eq_file(eqFile, sample)
    return ecCollection

## Main
def main():
    # User input
    usage = """%(prog)s is ... WIP.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-e", dest="equivalenceClassFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more salmon eq_classes.txt files")
    p.add_argument("-s", dest="sampleNames",
                   nargs="+",
                   required=True,
                   help="Input one or more sample names (paired to the equivalence classes)")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    
    args = p.parse_args()
    validate_args(args)
    
    ecCollection = parse_equivalence_classes(args.equivalenceClassFiles, args.sampleNames)
    
    ## TBD: Generate output from here
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
