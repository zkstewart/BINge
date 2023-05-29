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
                       the transcript ID you want to get counts for.
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        if isinstance(ecValue, int):
            return self.ec[ecValue]
        else:
            return self.ec[self.ids[ecValue]]
    
    def __repr__(self):
        return "<EquivalenceClassCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

class QuantCollection():
    def __init__(self):
        self.samples = [] # sample names
        self.quant = {} # quant count dicts
        
        self.numTranscripts = None # for checking quant file compatibility
    
    def parse_quant_file(self, quantFile, sample):
        # Validate input parameters
        assert os.path.isfile(quantFile), \
            f"Cannot parse '{quantFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Parse the file proper
        firstLine = True
        thisQuant = {}
        with open(quantFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                
                # Handle header line
                if firstLine == True:
                    assert sl == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"], \
                        "Quant file appears to lack the expected header? Cannot parse."
                    firstLine = False
                # Handle content lines
                else:
                    name, _, _, _, numReads = sl # might error here if file format is bad
                    thisQuant[name] = float(numReads)
        
        # Check compatibility of quant file with existing ones
        if self.numTranscripts == None:
            self.numTranscripts = len(thisQuant)
        else:
            assert self.numTranscripts == len(thisQuant), \
                "Quant files are incompatible since transcript numbers differ!"
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store counts inside parent structure
        for name, numReads in thisQuant.items():
            self.quant.setdefault(name, [])
            self.quant[name].append(float(numReads))
    
    def get_transcript_count(self, seqID):
        '''
        Parameters:
            seqID -- a string identifying the transcript ID you want to get counts for.
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        return self.quant[seqID]
    
    def __repr__(self):
        return "<QuantCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

# Define functions
def validate_salmon_files(salmonFiles):
    '''
    Validates Salmon files for 1) their existence and 2) their consistency of file format.
    Quits program if validation fails, so be warned!
    
    Parameters:
        salmonFiles -- a list containing strings pointing to Salmon files.
    Returns:
        fileFormat -- a string equal to "ec" if input files are equivalence classes,
                      or "quant" if they are quant.sf files.
    '''
    # Validate that salmon files exist
    for salmonFile in salmonFiles:
        if not os.path.isfile(salmonFile):
            print(f'I am unable to locate the salmon input file ({salmonFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate that input files are all of a consistent format
    isEC = False
    isQuant = False
    for salmonFile in salmonFiles:
        thisFileValid = False
        with open(salmonFile, "r") as fileIn:
            # Get the first 3 lines out of the file
            firstLine = fileIn.readline().rstrip("\r\n ")
            secondLine = fileIn.readline().rstrip("\r\n ")
            thirdLine = fileIn.readline().rstrip("\r\n ")
            
            # Check if it conforms to equivalence class expectations
            if firstLine.isdigit() and secondLine.isdigit() and not thirdLine.isdigit():
                isEC = True
                thisFileValid = True
            
            # Check if it conforms to quant file expectations
            elif firstLine.split("\t") == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
                isQuant = True
                thisFileValid = True
        
        if not thisFileValid:
            print(f"The input file '{salmonFile}' does not appear to be a Salmon quant or " + 
                  "equivalence class file")
            print("You should check your inputs and try again.")
            quit()
    
    if isEC and isQuant:
        print("You appear to have given a mix of quant and equivalence class files.")
        print("That's too hard for me to figure out, so please only give one type and try again.")
        quit()
    
    return "ec" if isEC else "quant"

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.bingeResultFile):
        print(f'I am unable to locate the BINge result file ({args.bingeResultFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate salmon files and set internal flag to mark what type of file we're using
    args.fileFormat = validate_salmon_files(args.salmonFiles)
    
    # Validate the input files are logically sound
    if len(args.salmonFiles) != len(args.sampleNames):
        print("Your salmon files and sample names are incompatible!")
        print((f"I'm seeing {len(args.salmonFiles)} salmon files " +
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

def parse_quants(quantFiles, sampleNames):
    '''
    Parses in one or more quant files from Salmon, producing a QuantCollection
    object enabling read count summarisation.
    
    Parameters:
        quantFiles -- a list containing strings pointing to the location
                      of Salmon quant files (quant.sf files).
        sampleNames -- an equal length list indicating the sample names for each
                       salmon quant file.
    '''
    assert len(quantFiles) == len(sampleNames), \
        ("parse_quants cannot parse quant files since the number " +
         f"of files ({len(quantFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    quantCollection = QuantCollection()
    for i in range(len(quantFiles)):
        quantFile = quantFiles[i]
        sample = sampleNames[i]
        
        quantCollection.parse_quant_file(quantFile, sample)
    return quantCollection

def parse_binge_clusters(bingeFile, typeToReturn="all"):
    '''
    Reads in the output file of BINge as a dictionary assocating clusters to their
    sequence members.
    
    Parameters:
        bingeFile -- a string pointing to the location of a BINge cluster output file.
        typeToReturn -- a string indicating whether to return all clusters ("all"), only
                        the binned clusters ("binned"), or only the unbinned clusters
                        ("unbinned")
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                             0: [seqid1, seqid2, ...],
                             1: [ ... ],
                             ...
                         }
    '''
    assert typeToReturn in ["all", "binned", "unbinned"]
    
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
                assert sl == ["cluster_num", "sequence_id", "cluster_type"], \
                    ("BINge file is expected to have a specific header line on the second line! " +
                     "Your file is hence not recognised as a valid BINge cluster file.")
                lineNum += 1
            
            # Handle content lines
            else:
                clustNum, seqID, clusterType = int(sl[0]), sl[1], sl[2]
                if typeToReturn == "all" or typeToReturn == clusterType:
                    clusterDict.setdefault(clustNum, [])
                    clusterDict[clustNum].append(seqID)
    return clusterDict

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside either 1) the equivalence class output files
    from Salmon, or 2) the quant.sf files output by Salmon. From either of these files, it
    will generate a tabular output of read counts per cluster suitable for DGE analysis.
    
    For each salmon file (-s), you should provide the sample name (-n). The order of these
    values should be equivalent, and will be reflected in the header of the output TSV file.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-s", dest="salmonFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more salmon eq_classes.txt or quant.sf files")
    p.add_argument("-n", dest="sampleNames",
                   nargs="+",
                   required=True,
                   help="Input one or more sample names (paired to the equivalence classes)")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    # Optional
    p.add_argument("--only_binned", dest="onlyBinned",
                   required=False,
                   action="store_true",
                   help="Optionally, only output counts for binned sequence clusters",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeResultFile,
        "all" if args.onlyBinned == False else "binned"
    )
    
    # Parse all salmon files
    if args.fileFormat == "ec":
        salmonCollection = parse_equivalence_classes(args.salmonFiles, args.sampleNames)
    else:
        salmonCollection = parse_quants(args.salmonFiles, args.sampleNames)
    
    # Generate output for clusters
    misses = 0
    numSeqs = 0
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("\t{0}\n".format("\t".join(args.sampleNames)))
        
        # Write cluster count lines
        for clusterNum, clusterIDs in clusterDict.items():
            
            # Sum counts across transcripts in this cluster per-sample
            clusterCount = [ 0 for _ in range(len(args.sampleNames)) ]
            for seqID in clusterIDs:
                numSeqs += 1
                try:
                    counts = salmonCollection.get_transcript_count(seqID)
                    clusterCount = [ clusterCount[x] + counts[x] for x in range(len(counts)) ]
                except: # this happens if Salmon filtered something out
                    misses += 1
                    continue
                
            # Output count
            fileOut.write("{0}\t{1}\n".format(
                f"cluster-{clusterNum}",
                "\t".join(map(str, clusterCount))
            ))
    
    # Alert user to potential problem
    if misses > 0:
        print("Potential non-issue warning:")
        print(" > Salmon filters out duplicate transcripts internally.")
        print(" > As a consequence, this program does expect there to be some differences between " +
              "the BINge clusters file and the equivalence class / quant files.")
        print(f" > In this case, it looks like Salmon filtered out {misses} sequences.")
        print(f" > Your BINge file indicates that {numSeqs} sequences were in your transcriptome.")
        print(" > If the number of filtered sequences is very similar to your total number of " +
              "sequences there might be a problem with your files.")
        print(" > i.e., your files were incompatible?")
        print(" > Otherwise no worries!")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
