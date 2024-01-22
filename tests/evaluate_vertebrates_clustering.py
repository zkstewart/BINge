#! python3
# evaluate_vertebrates_clustering.py

# This script is an adaptation of evaluate_clustering.py specifically for a
# BINge paper analysis. It is intended to facilitate a clustering evaluation
# where the true clusters are derived from the human-mouse-rat-zebrafish
# homolog groups.

import os, argparse, sys, re
import networkx as nx
from sklearn.metrics.cluster import adjusted_rand_score, rand_score, \
    normalized_mutual_info_score, adjusted_mutual_info_score

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.validation import validate_cluster_file
from modules.parsing import parse_binge_clusters

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.geneOrthologsFile):
        print(f'I am unable to locate the orthologs file ({args.geneOrthologsFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.clusterFile):
        print(f'I am unable to locate the cluster file ({args.clusterFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Make sure FASTA or GFF3 file exists
    for gff3File in args.annotationGFF3:
        if not os.path.isfile(gff3File):
            print(f'I am unable to locate the GFF3 file ({gff3File})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate input file format
    isBinge = validate_cluster_file(args.clusterFile)
    if not isBinge:
        raise ValueError("The input file was not validated as a BINge file!")
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_gff3_geneids(refseqGFF3Files):
    '''
    Parses a RefSeq annotations GFF3 from NCBI to a dictionary structure.
    
    Parameters:
        refseqGFF3Files -- a list containing strings that indicate the
                           location of  GFF3 files(s).
    Returns:
        idMappingDict -- a dictionary with structure like:
                    {
                        'seqID1': 'entrezID1'
                        'seqID2': 'entrezID2',
                        ...
                    }
    '''
    geneIDRegex = re.compile(r"GeneID:(\d+)")
    
    idMappingDict = {}
    for refseqGFF3File in refseqGFF3Files:
        with open(refseqGFF3File, "r") as fileIn:
            for line in fileIn:
                if not line.startswith("#"):
                    sl = line.rstrip("\r\n ").split("\t")
                    annotType, attributes = sl[2], sl[8]
                    if annotType == "mRNA":
                        seqID = attributes.split("ID=")[1].split(";")[0]
                        geneID = geneIDRegex.search(attributes).group(1)
                        
                        idMappingDict[seqID] = geneID
    return idMappingDict

def replace_test_with_entrez(testDict, idMappingDict):
    # Count how many occurrences there are of each Entrez ID
    countDict = {}
    for seqIDs in testDict.values():
        for seqID in seqIDs:
            entrezID = idMappingDict[seqID]
            countDict.setdefault(entrezID, 0)
            countDict[entrezID] += 1
    
    # Update any entrez IDs that are duplicated
    fixedIDs = {}
    for clustNum, seqIDs in testDict.items():
        newClust = []
        for seqID in seqIDs:
            # Get the entrez ID for this sequence
            entrezID = idMappingDict[seqID]
            
            # Fix the ID if it is duplicated
            if countDict[entrezID] > 1:
                fixedIDs.setdefault(entrezID, 0)
                fixedIDs[entrezID] += 1
                
                newClust.append(f"{entrezID}_{fixedIDs[entrezID]}")
            else:
                newClust.append(entrezID)
        
        # Store the fixed clusters
        testDict[clustNum] = newClust
    
    return testDict, fixedIDs

def parse_orthologs(orthologsFile, fixedIDs):
    '''
    Parses a gene_orthologs file from NCBI to a dictionary structure.
    
    Parameters:
        orthologsFile -- a string indicating the location of the gene_orthologs file.
        fixedIDs -- a dictionary containing entrez IDs which need to be substituted with
                    the amount of _${num} as in the value pairing.
    Returns:
        trueDict -- a dictionary with structure like:
                    {
                        'seqID1': 0, # number represents cluster ID
                        'seqID2': 1,
                        ...
                    }
        numGeneClusters -- an integer representing the number of clusters in the file.
    '''
    trueDict = {}
    numGeneClusters = -1
    
    with open(orthologsFile, "r") as fileIn:
        # Get first line to check file format
        firstLine = next(fileIn).rstrip("\r\n ")
        
        # Handle gene_orthologs file
        if firstLine.startswith("#tax_id\tGeneID"):
            graph = nx.Graph()
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                geneID, otherGeneID = sl[1], sl[4]
                
                geneIDs = [geneID] if not geneID in fixedIDs \
                    else [ f"{geneID}_{i+1}" for i in range(fixedIDs[geneID]) ]
                otherGeneIDs = [otherGeneID] if not otherGeneID in fixedIDs \
                    else [ f"{otherGeneID}_{i+1}" for i in range(fixedIDs[otherGeneID]) ]
                
                for gID in geneIDs:
                    graph.add_node(gID)
                for ogID in otherGeneIDs:
                    graph.add_node(ogID)
                
                for gID in geneIDs:
                    for ogID in otherGeneIDs:
                        graph.add_edge(gID, ogID)
            
            for connectedSeqIDs in nx.connected_components(graph):
                numGeneClusters += 1
                for seqID in connectedSeqIDs:
                    trueDict[seqID] = numGeneClusters
        
        # Raise error for unhandled file
        else:
            raise ValueError("The input file does not appear to be a gene_orthologs file!")
    
    return trueDict, numGeneClusters

## Main
def main():
    # User input
    usage = """%(prog)s is a modified evaluate_clustering.py script designed to receive not a
    reference annotation file, but the gene_orthologs file from https://ftp.ncbi.nih.gov/gene/DATA.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-g", dest="geneOrthologsFile",
                   required=True,
                   help="Input gene_orthologs file")
    p.add_argument("-c", dest="clusterFile",
                   required=True,
                   help="Input the CD-HIT or BINge cluster file")
    p.add_argument("-a", dest="annotationGFF3",
                   required=True,
                   nargs="+",
                   help="Input the RefSeq GFF3 file(s)")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for text results")
    # Optional
    p.add_argument("--beTolerant", dest="beTolerant",
                   required=False,
                   action="store_true",
                   help="""Optionally, if you are finding that this script errors out when 
                   comparing the number of sequences in your GFF3 to your cluster file, you
                   can provide this flag to prevent the error and make this program tolerant
                   to differences between the GFF3 and cluster file. You should only do this
                   if you find that the errors are not significant!""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the BINge cluster file
    testDict = parse_binge_clusters(args.clusterFile)
    
    # Parse the RefSeq GFF3 files
    idMappingDict = parse_gff3_geneids(args.annotationGFF3)
    
    # Drop any sequences from testDict that aren't in our idMappingDict
    """This can enable us to perform pairwise comparison of species by
    choosing which files to include with annotationGFF3"""
    keysToDelete = []
    for clustNum, seqIDs in testDict.items():
        testDict[clustNum] = [seqID for seqID in seqIDs if seqID in idMappingDict]
        if testDict[clustNum] == []:
            keysToDelete.append(clustNum)
    for key in keysToDelete:
        del testDict[key]
    
    # Replace test dict IDs with Entrez IDs, handling duplicates
    testDict, fixedIDs = replace_test_with_entrez(testDict, idMappingDict)
    
    # Parse the orthologs file
    trueDict, numGeneClusters = parse_orthologs(args.geneOrthologsFile, fixedIDs)
    
    # Flip the dict around and +numGeneClusters to prevent cluster number overlap
    testDict = {
        seqID : clustNum+numGeneClusters+1
        for clustNum, idList in testDict.items()
        for seqID in idList
    }
    
    # Drop any sequences in testDict that aren't in our trueDict
    "Must have the exact same sequences for comparison"
    for seqID in list(testDict.keys()):
        if not seqID in trueDict:
            del testDict[seqID]
    
    # Ensure that things are okay, erroring out if they are not
    if len(testDict) != len(trueDict):
        print("I wasn't able to make the cluster and GFF3 files match each other.")
        print("Specifically, I see {0} testable mRNA features in the GFF3".format(len(trueDict)))
        print(f"I see {len(testDict)} sequences in the cluster file (after pruning to match the GFF3)")
        # Change behaviour depending on whether we're tolerant or not
        if args.beTolerant:
            # Drop any sequences in trueDict that aren't in our testDict
            for seqID in list(trueDict.keys()):
                if not seqID in testDict:
                    del trueDict[seqID]
            
            print("--beTolerant was specified, which means I've tried to make these numbers match.")
            print(f"I've subset testable mRNAs in the GFF3 down to {len(trueDict)} features")
            
            if len(trueDict) == 0:
                print("After doing this, I now find zero (0) sequences to test.")
                print("Either the reference GFF3 is unusual, or your cluster file doesn't match?")
                print("Unable to continue from this point, so program will exit now.")
                quit()
            elif len(testDict) == len(trueDict):
                print("These numbers now match, so I can continue the evaluation.")
                print("But, if these numbers look odd, then something may have gone wrong.")
            else:
                print("These numbers still do not match, so I cannot continue with the evaluation.")
                print("This might be a problem with your input files, or maybe I have a bug...")
                print("If the former, fix things and try again. If the latter - sorry.")
                quit()
        else:
            print("--beTolerant was not specficied; hence, unless these numbers are identical, " +
                  "an evaluation cannot occur.")
            print("This might be a problem with your input files, or maybe I have a bug...")
            print("If the former, fix things and try again. If the latter - sorry.")
            quit()
    elif len(testDict) == 0:
        print("I'm not sure why, but I've found zero (0) sequences to test.")
        print("If this helps to debug, I see {0} testable mRNA features in the GFF3".format(len(trueDict)))
        print("Either the reference GFF3 is unusual, or your cluster file doesn't match?")
        print("Unable to continue from this point, so program will exit now.")
        quit()
    
    # Derive clustering labels for comparison
    trueList = []
    testList = []
    for mrnaID in trueDict.keys():
        trueList.append(trueDict[mrnaID])
        testList.append(testDict[mrnaID])
    
    # Score the clustering result
    arscore = adjusted_rand_score(trueList, testList)
    rscore = rand_score(trueList, testList)
    nmiscore = normalized_mutual_info_score(trueList, testList)
    amiscore = adjusted_mutual_info_score(trueList, testList)
    
    # Write and print output
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#evaluate_vertebrates_clustering.py clustering comparison results\n")
        
        # Write details
        gff3Files = ",".join(args.annotationGFF3)
        
        fileOut.write(f"Orthologs file\t{args.geneOrthologsFile}\n")
        fileOut.write(f"GFF3 file(s)\t{gff3Files}\n")
        fileOut.write(f"Cluster file\t{args.clusterFile}\n")
        fileOut.write(f"Number of testable mRNAs in annotation\t{len(trueList)}\n")
        fileOut.write(f"Number of testable genes in annotation\t{numGeneClusters}\n")
        fileOut.write(f"Number of predicted clusters\t{len(set(testList))}\n")
        fileOut.write(f"Rand Index Score\t{rscore}\n")
        fileOut.write(f"Adjusted Rand Index Score\t{arscore}\n")
        fileOut.write(f"NMI Score\t{nmiscore}\n")
        fileOut.write(f"AMI Score\t{amiscore}\n")
    
    print(f"Rand Index Score = {rscore}; Adjusted Rand Index Score = {arscore}; " + 
          f"NMI score = {nmiscore}; AMI score = {amiscore}; see output file for more details")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
