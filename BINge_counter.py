#! python3
# BINge_counter.py
# BIN Genes for Expression analyses - counting module

# Utility program for counting reads from Salmon equivalence classes,
# producing a tabulated output suitable for DGE analysis.

import os, argparse

from modules.validation import validate_salmon_files
from modules.parsing import parse_equivalence_classes, parse_quants, parse_binge_clusters

# Define functions
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
