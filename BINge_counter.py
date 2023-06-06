#! python3
# BINge_counter.py
# BIN Genes for Expression analyses - counting module

# Utility program for counting reads from Salmon equivalence classes,
# producing a tabulated output suitable for DGE analysis.

import os, argparse
from pathlib import Path

from modules.validation import validate_salmon_files
from modules.parsing import parse_dge_quants, parse_binge_clusters

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
    args.outputFileNames = []
    for suffix in [".length", ".counts", ".abundance"]:
        outputFileName = args.outputPrefix.rstrip("._ ") + suffix
        if os.path.isfile(outputFileName):
            print(f'File already exists at output location ({outputFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()
        args.outputFileNames.append(outputFileName)

def present_r_instructions(lengthFile, countFile, abundanceFile):
    '''
    Prints ( / presents) instructions for how to load in the resulting files
    from this program into R.
    
    Parameters:
        lengthFile -- a string indicating the location of the .length file
        countFile -- a string indicating the location of the .counts file
        abundanceFile -- a string indicating the location of the .abundance file
    '''
    # Convert all paths to posix format
    lengthFile = Path(os.path.abspath(lengthFile)).as_posix()
    countFile = Path(os.path.abspath(countFile)).as_posix()
    abundanceFile = Path(os.path.abspath(abundanceFile)).as_posix()
    
    # Format and print commands needed to load data into R
    instructionFormat = f'''# Locate gene abundance files
countsFile <- '{countFile}'
lengthFile <- '{lengthFile}'
abundanceFile <- '{abundanceFile}'
\n\
# Create list for tximport
bingeTx = list(
    counts = as.matrix(read.table(file=countsFile, header = TRUE, sep = "\\t", stringsAsFactors = FALSE, row.names = 1)),
    abundance = as.matrix(read.table(file=abundanceFile, header = TRUE, sep = "\\t", stringsAsFactors = FALSE, row.names = 1)),
    length = as.matrix(read.table(file=lengthFile, header = TRUE, sep = "\\t", stringsAsFactors = FALSE, row.names = 1)),
    countsFromAbundance = "no"
)
\n\
# Create DESeq2 dataset via tximport
dds <- DESeqDataSetFromTximport(bingeTx, colData = coldata.table, design = ~ 1)'''
    
    return instructionFormat

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for use downstream of BINge clustering. It
    will read in the output of BINge alongside the quant.sf files output by Salmon;
    equivalence class files are not supported, since the EffectiveLength value in the
    quant.sf output is necessary for downstream DGE analysis. This script will produce
    three files ($PREFIX.length, $PREFIX.counts, $PREFIX.abundance) which can be loaded
    into R for use with DESeq2's DESeqDataSetFromTximport.
    
    For each salmon file (-s), you should provide the sample name (-n). The order of these
    values should be equivalent, and will be reflected in the header of the output files.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeResultFile",
                   required=True,
                   help="Input the TSV result of running BINge")
    p.add_argument("-s", dest="salmonFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more Salmon quant.sf files")
    p.add_argument("-n", dest="sampleNames",
                   nargs="+",
                   required=True,
                   help="Input one or more sample names (paired to the Salmon files)")
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Output file prefix for tximport-formatted results")
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
        print("BINge_counter does not support equivalence class input, sorry.")
        print("Program will exit now")
        quit()
    else:
        salmonDGECollection = parse_dge_quants(args.salmonFiles, args.sampleNames)
    
    # Generate output for clusters
    misses = 0
    numSeqs = 0
    zeros = 0
    lengthFile, countFile, abundanceFile = args.outputFileNames
    
    with open(lengthFile, "w") as lengthOut, open(countFile, "w") as countOut, open(abundanceFile, "w") as abundanceOut:
        # Write headers
        lengthOut.write("\t{0}\n".format("\t".join(args.sampleNames)))
        countOut.write("\t{0}\n".format("\t".join(args.sampleNames)))
        abundanceOut.write("\t{0}\n".format("\t".join(args.sampleNames)))
        
        # Write cluster lines
        for clusterNum, clusterIDs in clusterDict.items():
            
            # Sum counts, abundance, and length across transcripts in this cluster per-sample
            clusterCount = [ 0 for _ in range(len(args.sampleNames)) ]
            clusterAbundance = [ 0 for _ in range(len(args.sampleNames)) ]
            for seqID in clusterIDs:
                numSeqs += 1
                try:
                    counts = salmonDGECollection.get_transcript_count(seqID)
                    clusterCount = [ clusterCount[x] + counts[x] for x in range(len(counts)) ]
                    
                    abundances = salmonDGECollection.get_transcript_tpm(seqID)
                    clusterAbundance = [ clusterAbundance[x] + abundances[x] for x in range(len(abundances)) ]
                except: # this happens if Salmon filtered something out
                    misses += 1
                    continue
            
            # Calculate weighted average of effective transcript length
            "This tries to replicate https://rdrr.io/bioc/tximport/src/R/summarizeToGene.R"
            if sum(clusterAbundance) > 0:
                # First, for this "row", get the weighted length per transcript
                "rowsum would do the below across an entire dataframe at once; we do it per-row"
                weightedLengths = [ 0 for _ in range(len(args.sampleNames)) ]
                for seqID in clusterIDs:
                    try:
                        abundances = salmonDGECollection.get_transcript_tpm(seqID)                        
                        lengths = salmonDGECollection.get_transcript_effective_length(seqID)
                        wl = [ abundance * length for abundance, length in zip(abundances, lengths) ]
                        weightedLengths = [ weightedLengths[x] + wl[x] for x in range(len(wl)) ]
                    except:
                        continue
                
                # Get the length offset [borrowing lengthMat variable name from R source code]
                lengthMat = [ wl / abundance if not (wl == 0 or abundance == 0) else None for wl, abundance in zip(weightedLengths, clusterAbundance) ]
                
                # Impute missing values
                "DESeq2 uses the geometric mean with tximport, so we'll do the same"
                geometricMean = 1
                numNotMissing = 0
                for l in lengthMat:
                    if l != None:
                        geometricMean *= l
                        numNotMissing += 1
                geometricMean = geometricMean**(1/numNotMissing)
                
                lengthMat = [ l if l != None else geometricMean for l in lengthMat ]
            
            # If this cluster has no read depth, exclude it from output
            else:
                zeros += 1
                continue
            
            # Output to each file
            lengthOut.write("{0}\t{1}\n".format(
                f"cluster-{clusterNum}",
                "\t".join(map(str, lengthMat))
            ))
            countOut.write("{0}\t{1}\n".format(
                f"cluster-{clusterNum}",
                "\t".join(map(str, clusterCount))
            ))
            abundanceOut.write("{0}\t{1}\n".format(
                f"cluster-{clusterNum}",
                "\t".join(map(str, clusterAbundance))
            ))
    
    # Provide basic statistics to user which can be helpful for identifying issues
    print("# BINge_counter statistics on input clusters")
    print(f" > Number of clusters = {len(clusterDict)}.")
    print(f" > Number of transcripts = {numSeqs}.\n")
    
    # Alert user to potential problem(s)
    if misses > 0 or zeros > 0:
        print("# Potential problem warning(s):")
        if misses > 0:
            print(" > MISSES:")
            print(" > Salmon filters out duplicate transcripts internally.")
            print(" > As a consequence, this program does expect there to be some differences between " +
                "the BINge clusters file and the equivalence class / quant files.")
            print(f" > In this case, it looks like Salmon filtered out {misses} sequences.")
            print(" > If the number of filtered sequences is very similar to your total number of " +
                "sequences there might be a problem with your files.")
            print(" > i.e., your files were incompatible?")
            print(" > Otherwise no worries!\n")
        if zeros > 0:
            print(" > ZEROS:")
            print(" > Some clusters in your data are associated with 0 read depth.")
            print(f" > In this case, it looks like {zeros} clusters had no read alignment.")
            print(" > FYI - these clusters were NOT written to file.")
            print(" > If the number of clusters with no read alignment is very similar to " +
                "your total number of clusters, there might be a problem with your files.")
            print(" > i.e., your read alignments are of poor quality?")
            print(" > Otherwise no worries!\n")
    
    # Provide instruction for loading data into R
    print("# BINge_counter: how to use these outputs in R")
    print(" > Adapt the following code to load these outputs into R")
    print(" > Load in the coldata and set the design as appropriate")
    print("\n\n##########")
    print(present_r_instructions(lengthFile, countFile, abundanceFile))
    print("##########\n\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
