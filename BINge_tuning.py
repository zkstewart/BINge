#! python3
# BINge_tuning.py
# BIN Genes for Expression analyses - CD-HIT tuning module

# Utility program for tuning CD-HIT parameters prior to clustering of
# unbinned sequences.

import os, argparse, sys
from sklearn.metrics.cluster import adjusted_rand_score

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts import ZS_GFF3IO, ZS_ClustIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print(f'I am unable to locate the GFF3 annotation file ({args.gff3File})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.clstrFile):
        print(f'I am unable to locate the CD-HIT .clstr file ({args.clstrFile})')
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
    usage = """%(prog)s is a module intended for use prior to BINge clustering. The final
    step of BINge clustering is to use CD-HIT to cluster any unbinned sequences de novo,
    without using the genome. Having CD-HIT correctly tuned to render good results is
    important for the reliability of its results.
    
    This module allows one to test CD-HIT's clustering on an existing genome annotation.
    Provide a GFF3 file for a species related to yours which has a high-quality annotation
    inclusive of alternatively spliced isoforms. Additionally, provide the .clstr output
    of CD-HIT that you've run beforehand on the isoform transcripts for that same genome.
    This script will then assess how well CD-HIT was able to re-discover the gene groupings.
    
    Note that the identity value will be difficult to assess since multi-subspecies use will
    necessitate a lower identity threshold than clustering same species genes.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input the GFF3 annotation file")
    p.add_argument("-c", dest="clstrFile",
                   required=True,
                   help="Input the CD-HIT output .clstr file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for text results")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the CD-HIT cluster file
    testDict = ZS_ClustIO.CDHIT.parse_clstr_file(args.clstrFile) # we will 'test' this against our ground truth
    
    # Parse the GFF3 file into memory
    gff3 = ZS_GFF3IO.GFF3(args.gff3File)
    
    # Derive a cluster dict from the GFF3 to use as our true labels
    clusterNum = 1
    trueDict = {} # stores our 'true' cluster assignments based on the annotation itself
    for geneFeature in gff3.types["gene"]:
        if hasattr(geneFeature, "mRNA"):
            for mrnaFeature in geneFeature.mRNA:
                trueDict[mrnaFeature.ID] = clusterNum
            clusterNum += 1
    
    # Drop any sequences in testDict that aren't in our trueDict
    "Must have the exact same sequences for comparison"
    foundNum = 0
    for seqID in list(testDict.keys()):
        if not seqID in trueDict:
            del testDict[seqID]
        foundNum += 1
    
    # Ensure that things are okay, erroring out if they are not
    if len(testDict) != len(gff3.types["mRNA"]):
        print("I wasn't able to make the cluster and GFF3 files match each other.")
        print("Specifically, I see {0} mRNA features in the GFF3".format(len(gff3.types["mRNA"])))
        print(f"I see {len(testDict)} sequences in the cluster file (after pruning to match the GFF3)")
        print("Unless these numbers are identical, an evaluation cannot occur.")
        print("This might be a problem with your input files, or maybe I have a bug...")
        print("If the former, fix things and try again. If the latter - sorry.")
        quit()
    
    # Derive clustering labels for comparison
    trueList = []
    testList = []
    for mrnaFeature in gff3.types["mRNA"]:
        trueList.append(trueDict[mrnaFeature.ID])
        testList.append(testDict[mrnaFeature.ID])
    
    # Score the clustering result
    score = adjusted_rand_score(trueList, testList)
    
    # Write and print output
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge_tuning clustering comparison\n")
        
        # Write details
        fileOut.write(f"Annotation file\t{args.gff3File}\n")
        fileOut.write(f"Cluster file\t{args.clstrFile}\n")
        fileOut.write(f"Number of sequences\t{len(trueList)}\n")
        fileOut.write(f"Number of genes in annotation\t{len(set(trueList))}\n")
        fileOut.write(f"Number of clusters from CD-HIT\t{len(set(testList))}\n")
        fileOut.write(f"Adjusted Rand Index Score\t{score}\n")
    
    print(f"Adjusted Rand Index Score = {score}; see output file for more details")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
