#! python3
# cdhit_parameter_pipeline.py

# Utility program for determining the optimal combination of
# CD-HIT parameters to recapitulate a reference genome's annotation
# by using BINge_tuning.py. Specifically, this script is run BEFORE
# BINge_tuning.py, and is intended to automate the process of trialing
# many parameter combinations for testing to see which is most optimal
# for the given reference genome.

import os, argparse, sys
from itertools import product
from threading import Thread

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.cdhit import CDHIT

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.mrnaFastaFile):
        print(f'I am unable to locate the input FASTA file ({args.mrnaFastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.cdhitDir):
        print(f'I am unable to locate the CD-HIT executable directory ({args.cdhitDir})')
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    if not (os.path.isfile(os.path.join(args.cdhitDir, "cd-hit")) or os.path.isfile(os.path.join(args.cdhitDir, "cd-hit.exe"))):
        print(f"I am unable to locate the cd-hit executable at '{args.cdhitDir}'")
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    if not (os.path.isfile(os.path.join(args.cdhitDir, "cd-hit-est")) or os.path.isfile(os.path.join(args.cdhitDir, "cd-hit-est.exe"))):
        print(f"I am unable to locate the cd-hit-est executable '{args.cdhitDir}'")
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    # Validate CD-HIT parameters
    for x in args.identity:
        if args.isProtein:
            if not 0.40 <= x <= 1.0:
                print(f"For proteins, identity value must be in the range 0.40 <= {x} <= 1.0")
                quit()
        else:
            if not 0.80 <= x <= 1.0:
                print(f"For nucleotides, identity value must be in the range 0.80 <= {x} <= 1.0")
                quit()
    for x in args.shortCov:
        if not 0 < x <= 1.0:
            print(f"shortCov value must be in the range 0 < {x} <= 1.0")
            quit()
    for x in args.longCov:
        if not 0 < x <= 1.0:
            print(f"longCov value must be in the range 0 < {x} <= 1.0")
            quit()
    if not 0 < args.paramThreads:
        print("paramThreads value must be an integer greater than 0")
        quit()
    if not 0 < args.cdhitThreads:
        print("cdhitThreads value must be an integer greater than 0")
        quit()
    # Generate parameter combinations
    args.parameters = list(product(
        args.identity, args.shortCov, args.longCov, [args.useLocal]))
    
    # Generate output file names based on parameter combinations
    args.outputPrefixes = [
        os.path.join(args.outputDirectory, 
                     "_".join(param))
        for param in list(product([ f"c{identity}" for identity in args.identity ],
                         [ f"aS{shortCov}" for shortCov in args.shortCov ],
                         [ f"aL{longCov}" for longCov in args.longCov ],
                         [ "G0" if args.useLocal else "G1"]))
    ]
    args.outputFileNames = [
        prefix + ".clstr"
        for prefix in args.outputPrefixes
    ]
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory):
        print(f'Output location already exists ({args.outputDirectory})')
        for fileName in args.outputFileNames:
            if os.path.exists(fileName):
                print(f"'{fileName}', which would be created, already exists.")
                print("Set a different output directory, or move the existing files, then try again.")
                quit()
        print("It doesn't look like things will conflict, so I'll keep moving ahead.")
    elif os.path.exists(args.outputDirectory):
        print(f'Something other than a directory already exists at the output location ({args.outputDirectory})')
        print('Make sure you specify a new directory, or move whatever is currently here then try again.')
        quit()
    else:
        os.mkdir(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' created as part of argument validation.")

class CDHITWorkerThread(Thread):
    '''
    This provides a modified Thread which encapsulates the code needed
    to run CD-HIT.
    '''
    def __init__(self, inputFileName, outputFilePrefix, cdhitDir, mem,
                 identity, shortCov, longCov, isLocal, threads, molecule):
        Thread.__init__(self)
        
        self.mem = mem
        self.identity = identity
        self.shortCov = shortCov
        self.longCov = longCov
        self.local = isLocal
        self.threads = threads
        self.molecule = molecule
        
        self.inputFileName = inputFileName
        self.outputFilePrefix = outputFilePrefix
        self.cdhitDir = cdhitDir
        
        self.resultClstrFile = None
        self.exception = None
    
    def run_cdhit(self):
        # Configure clusterer object
        clusterer = CDHIT(self.inputFileName, self.molecule, self.cdhitDir)
        clusterer.identity = self.identity
        clusterer.threads = 1
        clusterer.mem = self.mem
        clusterer.clean = False
        clusterer.threads = self.threads
        clusterer.set_shorter_cov_pct(self.shortCov)
        clusterer.set_longer_cov_pct(self.longCov)
        if self.local == True:
            clusterer.set_local()
        else:
            clusterer.set_global()
        
        # Get output details
        outputDir = os.path.dirname(os.path.abspath(self.outputFilePrefix))
        outputFasta = os.path.basename(os.path.abspath(self.outputFilePrefix))
        
        # Run clustering
        clusterer.cdhit(self.inputFileName, outputDir, outputFasta)
        
        # Clean up FASTA result but retain clstr file location
        os.unlink(os.path.abspath(self.outputFilePrefix))
        self.resultClstrFile = os.path.abspath(os.path.abspath(self.outputFilePrefix) + ".clstr")
    
    def run(self):
        try:
            self.run_cdhit()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception

## Main
def main():
    # User input
    usage = """%(prog)s will run CD-HIT many times using every possible parameter
    configuration as specified in this script. The result will be a .clstr output file
    per parameter configuration, which can be tested with BINGe_tuning.py to assess which
    parameters give the best clustering result.
    
    Each parameter (e.g., -c) should be given one or more values that you'd like to test.
    Every possible configuration of parameters will be tested, so the number of times
    CD-HIT will run scales very quickly.
    
    Note that the two thread arguments are multiplicative. If you specify --paramThreads 2,
    as well as --cdhitThreads 4, then this program will run two separate CD-HIT processes
    each with 4 threads; hence, there should be 8 cores available.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="mrnaFastaFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output file name")
    p.add_argument("-cdhit", dest="cdhitDir",
                   required=True,
                   help="""Specify the directory containing
                   the cd-hit and cd-hit-est executables""")
    p.add_argument("-c", dest="identity",
                   required=True,
                   type=float,
                   nargs="+",
                   help="""Specify one or more -c parameters to provide
                   CD-HIT""")
    p.add_argument("-aS", dest="shortCov",
                   required=True,
                   type=float,
                   nargs="+",
                   help="""Specify one or more -aS parameters to provide
                   CD-HIT""")
    p.add_argument("-aL", dest="longCov",
                   required=False,
                   type=float,
                   nargs="+",
                   help="""Specify one or more -aL parameters to provide
                   CD-HIT""")
    # Optional
    p.add_argument("--isProtein", dest="isProtein",
                   required=False,
                   action="store_true",
                   help="""Optionally specify this parameter if your input
                   file contains protein sequences, as opposed to nucleotides""",
                   default=False)
    p.add_argument("--local", dest="useLocal",
                   required=False,
                   action="store_true",
                   help="""Optionally specify this parameter to run CD-HIT
                   using local alignment rather than the default global""",
                   default=False)
    p.add_argument("--mem", dest="mem",
                   required=False,
                   type=int,
                   help="""Specify how many megabytes of memory to
                   provide CD-HIT (default==6000)""",
                   default=6000)
    p.add_argument("--paramThreads", dest="paramThreads",
                   required=False,
                   type=int,
                   help="""Specify how many separate CD-HIT processes should
                   occur simultaneously (default==1)""",
                   default=1)
    p.add_argument("--cdhitThreads", dest="cdhitThreads",
                   required=False,
                   type=int,
                   help="""Specify how many threads to run each CD-HIT process
                   with (default==1)""",
                   default=1)
    
    args = p.parse_args()
    validate_args(args)
    
    # Run CD-HIT processing threads
    for i in range(0, len(args.parameters), args.paramThreads): # only process n (threads) files at a time
        processing = []
        for x in range(args.paramThreads): # begin processing n CD-HIT runs
            if i+x < len(args.parameters): # parent loop may excess if n > the number of param combinations
                identity, shortCov, longCov, isLocal = args.parameters[i+x]
                outputFilePrefix = args.outputPrefixes[i+x]
                
                cdhitWorkerThread = CDHITWorkerThread(args.mrnaFastaFile, outputFilePrefix,
                                                      args.cdhitDir, args.mem, identity,
                                                      shortCov, longCov, isLocal, args.cdhitThreads,
                                                      "protein" if args.isProtein else "nucleotide")
                processing.append(cdhitWorkerThread)
                cdhitWorkerThread.start()
        
        # Gather results
        for cdhitWorkerThread in processing:
            # Wait for thread to end
            cdhitWorkerThread.join()
            
            # Print success of this processing thread
            print(f"Finished processing {cdhitWorkerThread.resultClstrFile}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
