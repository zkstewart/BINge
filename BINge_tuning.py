#! python3
# BINge_tuning.py
# BIN Genes for Expression analyses - CD-HIT tuning module

# Utility program for tuning CD-HIT parameters prior to clustering of
# unbinned sequences.

import os, argparse, sys
from sklearn.metrics.cluster import adjusted_rand_score

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts import ZS_GFF3IO, ZS_ClustIO

from modules.validation import validate_cluster_file
from modules.parsing import parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print(f'I am unable to locate the GFF3 annotation file ({args.gff3File})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.clusterFile):
        print(f'I am unable to locate the cluster file ({args.clusterFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate input file format
    args.isBinge = validate_cluster_file(args.clusterFile)
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

## Main
def main():
    # User input
    usage = """%(prog)s is a module intended for two primary purposes.
    
    First, it can be used to fine-tune CD-HIT's parameters prior to BINge's clustering.
    Specifically, the final step of BINge clustering is to use CD-HIT to cluster any
    unbinned sequences de novo, without using the genome. Having CD-HIT correctly tuned
    to render good results is important for the reliability of its results.
    
    Second, it can be used to test CD-HIT's or BINge's clustering on an existing genome
    annotation. Provide a GFF3 file for a species related to yours which has a high-quality
    annotation inclusive of alternatively spliced isoforms. Additionally, provide the .clstr
    of CD-HIT or the .tsv output of BINge that you've run beforehand on the isoform transcripts
    for that same genome. This script will then assess how well the clusterer was able to
    re-discover the gene groupings.
    
    Note that (for the first purpose) the identity value will be difficult to assess since
    multi-subspecies use will necessitate a lower identity threshold than clustering same
    species genes. But the remaining parameters should hold true.
    
    Lastly, this program will automatically detect whether the given cluster file comes
    from BINge or CD-HIT.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input the GFF3 annotation file")
    p.add_argument("-c", dest="clusterFile",
                   required=True,
                   help="Input the CD-HIT or BINge cluster file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for text results")
    # Optional
    p.add_argument("--seq_prefix", dest="seqPrefix",
                   required=False,
                   help="""Optionally, if your GFF3 has a prefix before all its sequence
                   IDs that is not seen in your FASTA file of transcripts, specify that here.
                   For example, with NCBI GenBank genomes, it's not uncommon for the transcript
                   FASTA file to have IDs like 'XM_009121514.3' but the GFF3 would index that
                   as 'rna-XM_009121514.3'. So you would specify 'rna-' here to address that.""",
                   default="")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the GFF3 file into memory
    gff3 = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    
    # Derive a cluster dict from the GFF3 to use as our true labels
    numGeneClusters = 0
    trueDict = {} # stores our 'true' cluster assignments based on the annotation itself
    for geneFeature in gff3.types["gene"]:
        if hasattr(geneFeature, "mRNA"):
            for mrnaFeature in geneFeature.mRNA:
                trueDict[mrnaFeature.ID] = numGeneClusters
            numGeneClusters += 1
    
    # Parse the CD-HIT / BINge cluster file, changing cluster IDs to not overlap
    if args.isBinge:
        testDict = parse_binge_clusters(args.clusterFile) # we will 'test' this against our ground truth
    else:
        testDict = ZS_ClustIO.CDHIT.parse_clstr_file(args.clusterFile) 
    
    # Flip the dict around and +numGeneClusters to prevent cluster number overlap
    testDict = {
        f"{args.seqPrefix}{seqID}" : clustNum+numGeneClusters
        for clustNum, idList in testDict.items()
        for seqID in idList
    }
    
    # Drop any sequences in testDict that aren't in our trueDict
    "Must have the exact same sequences for comparison"
    for seqID in list(testDict.keys()):
        if not seqID in trueDict:
            del testDict[seqID]
    
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
        fileOut.write(f"Cluster file\t{args.clusterFile}\n")
        fileOut.write(f"Number of sequences\t{len(trueList)}\n")
        fileOut.write(f"Number of genes in annotation\t{numGeneClusters}\n")
        fileOut.write(f"Number of predicted clusters\t{len(set(testList))}\n")
        fileOut.write(f"Adjusted Rand Index Score\t{score}\n")
    
    print(f"Adjusted Rand Index Score = {score}; see output file for more details")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
