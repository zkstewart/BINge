#! python3
# tabulate_salmon_qc.py

# Utility program to obtain output information from Salmon
# mapping and tabulate some QC statistics. Specifically, it
# helps to identify 1) what the read mapping rate was initially,
# and 2) how many of those read mappings made it through the
# clustering process.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.parsing import parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    stderrDir = os.path.dirname(os.path.abspath(args.salmonStderrPrefix))
    stderrPrefix = os.path.basename(args.salmonStderrPrefix)
    if not os.path.isdir(stderrDir):
        print(f'I am unable to locate the parent directory where Salmon stderr files should be ({args.stderrDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not any([ True if f.startswith(stderrPrefix) else False for f in os.listdir(stderrDir) ]):
        print(f'I am unable to locate any files with "{stderrPrefix}" prefix at ({args.stderrDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    if not os.path.isdir(args.salmonOutputsDir):
        print(f'I am unable to locate the parent directory where Salmon subdirectories should be ({args.salmonOutputsDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    
    args.stderrDir = stderrDir
    args.stderrPrefix = stderrPrefix

def parse_salmon_stderr(stderrFiles):
    '''
    Parameters:
        stderrFiles -- a list containing strings pointing to the locations of 
                       salmon stderr files.
    Returns:
        stderrDict -- a dictionary with structure like:
                      {
                          'filename1': [ mappingRate, 'outputDir' ],
                          'filename2': [ ... ],
                          ...
                      }
    '''
    stderrDict = {}
    for file in stderrFiles:
        stderrDict[os.path.basename(file)] = [None, None]
        
        with open(file, "r") as fileIn:
            for line in fileIn:
                # Get mapping rate
                if "[info] Mapping rate =" in line:
                    mappingRate = line.rstrip("\r\n ").split(" = ")[1]
                    stderrDict[os.path.basename(file)][0] = mappingRate
                # Get output directory
                elif "### [ output ] => " in line:
                    outDir = line.rstrip("\r\n ").split("{ ")[1].split(" }")[0]
                    stderrDict[os.path.basename(file)][1] = outDir
        
        assert all([ v != None for v in stderrDict[os.path.basename(file)] ]), \
            f"Parsing salmon output failed; {os.path.basename(file)} = {stderrDict[os.path.basename(file)]}"
    
    return stderrDict

## Main
def main():
    # User input
    usage = """%(prog)s locates the outputs of Salmon and tabulates some relevant QC statistics
    pertaining to the read mapping success of each sample. It anticipates there being stderr files
    generated via running Salmon, as well as quant.sf files as output by Salmon in subdirectories.
    
    With relation to BINge, it allows one to identify what proportion of read mappings are
    represented by the clustered output. Excessive amounts of reads not being factored in may
    indicate a poor clustering output. Otherwise, such statistics may be necessary to report in a
    publication.
    """
    p = argparse.ArgumentParser(description=usage)
    
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-s", dest="salmonStderrPrefix",
                   required=True,
                   help="""Specify the prefix which uniquely identifies Salmon stderr (.e) files;
                   you may want to specify the path leading up to the files as well""")
    p.add_argument("-d", dest="salmonOutputsDir",
                   required=True,
                   help="Specify location where Salmon has written output subfolders to")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the QC statistics table")
    # Optional
    p.add_argument("--only_binned", dest="onlyBinned",
                   required=False,
                   action="store_true",
                   help="Optionally, consider only binned sequence clusters",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all .o files
    stderrFiles = [ os.path.join(args.stderrDir, f) for f in os.listdir(args.stderrDir) if f.startswith(args.stderrPrefix) ]
    
    # Parse .o files
    stderrDict = parse_salmon_stderr(stderrFiles)
    
    # Validate that all output directories indicated by stderr are locateable
    for key, rateAndDir in stderrDict.items():
        mappingRate, outDir = rateAndDir
        fullOutDir = os.path.join(args.salmonOutputsDir, outDir)
        if not os.path.isdir(fullOutDir):
            try:
                fullOutDir = os.path.abspath(outDir)
                assert os.path.isdir(fullOutDir)
            except:
                print(f"Expected to find Salmon outputs at {fullOutDir}, but directory doesn't exist!?")
                print("Check that your -d parameter is correct and try again.")
                quit()
        
        quantFile = os.path.join(fullOutDir, "quant.sf")
        if not os.path.isfile(quantFile):
            print(f"Expected to find quant.sf file at {fullOutDir}, but file appears to be absent!?")
            print("Check that you actually ran Salmon in quant mode and try again.")
            quit()
        
        # Set quant file location since its been validated and we may have entered the try clause earlier
        stderrDict[key] = [mappingRate, quantFile, outDir] # outDir gives us a better sample name
    
    # Parse the BINge cluster file
    clusterDict = parse_binge_clusters(
        args.bingeResultFile,
        "all" if args.onlyBinned == False else "binned"
    )
    
    # Extract set of sequence IDs which are part of a BINge cluster
    clusterSet = set().union(*clusterDict.values())
    
    # Parse Salmon output files
    salmonCountsDict = {}
    for _, rateQuantOutdir in stderrDict.items():
        _, quantFile, outDir = rateQuantOutdir # deconstructs [mappingRate, quantFile, outDir]
        salmonCountsDict[outDir] = [0, 0] # [num in clusters, num not in clusters]
        
        with open(quantFile, "r") as fileIn:
            firstLine = True
            for line in fileIn:
                if firstLine == True:
                    firstLine = False
                    continue
                else:
                    seqName, length, effectiveLength, tpm, \
                        numReads = line.rstrip("\r\n ").split("\t")
                    
                    if seqName in clusterSet:
                        salmonCountsDict[outDir][0] += float(numReads)
                    else:
                        salmonCountsDict[outDir][1] += float(numReads)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("{0}\n".format("\t".join([
            "#salmon_dir", "mapping_rate", "cluster_proportion", "clustered_mapping_rate"
        ])))
        
        # Write content lines
        for outDir, counts in salmonCountsDict.items():
            # Get the mapping rate for this sample
            mappingRate = float([
                _mappingRate
                for _mappingRate, _, _outDir in stderrDict.values()
                if _outDir == outDir
            ][0][:-1]) # trim the '%' symbol and make it a float
            
            # Calculate the proportion of counts represented in a cluster
            clusterProportion = 1 - (counts[1] / counts[0])
            clusterMapRate = mappingRate * clusterProportion
            
            # Write to file
            fileOut.write("{0}\n".format("\t".join([
                outDir, str(mappingRate), str(clusterProportion * 100),
                str(clusterMapRate)
            ])))
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
