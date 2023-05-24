#! python3
# BINge_counter.py
# BIN Genes for Expression analyses - counting module

# Utility program for counting reads from Salmon equivalence classes,
# producing a tabulated output suitable for DGE analysis.

import os, argparse

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
            thisEC = { i:0.0 for i in range(0, numTranscripts) } # holds onto counts per transcript for this file
            for _ in range(numECs):
                sl = list(map(int, fileIn.readline().rstrip("\r\n ").split("\t")))
                numIDs, idIndices, numReads = sl[0], sl[1:-1], sl[-1]
                numReads = numReads / numIDs
                
                # Store relevant data
                for idIndex in idIndices:
                    thisEC[idIndex] += numReads
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store EC counts inside parent structure
        for ecIndex, ecCount in thisEC.items():
            self.ec.setdefault(ecIndex, [])
            self.ec[ecIndex].append(ecCount)
    
    def get_transcript_count(self, ecValue):
        '''
        Parameters:
            ecValue -- an integer identifying the equivalence class index
                       you want to get counts for, OR a string identifying
                       the transcript ID you want to get counts for
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples
        '''
        if isinstance(ecValue, int):
            return self.ec[ecValue]
        else:
            return self.ec[self.ids[ecValue]]
    
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
    '''
    Parses in one or more equivalence class files from Salmon, producing
    an EquivalenceClassCollection object enabling read count summarisation.
    
    Parameters:
        equivalenceClassFiles -- a list containing strings pointing to the location
                                 of Salmon equivalence class files
                                 (eq_classes.txt files).
        sampleNames -- an equal length list indicating the sample names for each
                       equivalence class file.
    '''
    assert len(equivalenceClassFiles) == len(sampleNames), \
        ("parse_equivalence_classes cannot parse equivalence classes since the number " +
         f"of files ({len(equivalenceClassFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    ecCollection = EquivalenceClassCollection()
    for i in range(len(equivalenceClassFiles)):
        eqFile = equivalenceClassFiles[i]
        sample = sampleNames[i]
        
        ecCollection.parse_eq_file(eqFile, sample)
    return ecCollection

def parse_binge_clusters(bingeFile):
    '''
    Reads in the output file of BINge as a dictionary assocating clusters to their
    sequence members.
    
    Parameters:
        bingeFile -- a string pointing to the location of a BINge cluster output file.
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                             0: [seqid1, seqid2, ...],
                             1: [ ... ],
                             ...
                         }
    '''
    clusterDict = {}
    lineNum = 0
    with open(bingeFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if lineNum == 0:
                assert line.startswith("#BINge clustering information file"), \
                    ("BINge file is expected to start with a specific comment line! " +
                    "Your file is hence not recognised as a valid BINge cluster file.")
                lineNum += 1
            elif lineNum == 1:
                assert sl == ["cluster_num", "sequence_id"], \
                    ("BINge file is expected to have a specific header line on the second line! " +
                     "Your file is hence not recognised as a valid BINge cluster file.")
                lineNum += 1
            
            # Handle content lines
            else:
                clustNum, seqID = int(sl[0]), sl[1]
                clusterDict.setdefault(clustNum, [])
                clusterDict[clustNum].append(seqID)
    return clusterDict
                
## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside one or more equivalence class files from
    Salmon, generating a tabular output of read counts per cluster suitable for DGE analysis.
    
    For each equivalence class file (-e), you should provide the name of the sample (-s).
    The order of these values should be equivalent, and will be reflected in the header
    of the output TSV file.
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
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(args.bingeResultFile)
    
    # Parse all equivalence classes
    ecCollection = parse_equivalence_classes(args.equivalenceClassFiles, args.sampleNames)
    
    # Generate output for clusters
    misses = []
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("\t{0}\n".format("\t".join(args.sampleNames)))
        
        # Write cluster count lines
        for clusterNum, clusterIDs in clusterDict.items():
            
            # Sum counts across transcripts in this cluster per-sample
            clusterCount = [ 0 for _ in range(len(args.sampleNames)) ]
            for seqID in clusterIDs:
                try:
                    counts = ecCollection.get_transcript_count(seqID)
                    clusterCount = [ clusterCount[x] + counts[x] for x in range(len(counts)) ]
                except: # this probably happens if Salmon filtered something out?
                    misses.append(seqID)
                    continue
                
            # Output count
            fileOut.write("{0}\t{1}\n".format(
                f"cluster-{clusterNum}",
                "\t".join(map(str, clusterCount))
            ))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
