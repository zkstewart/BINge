#! python3
# evaluate_vertebrates_clustering.py

# This script is an adaptation of evaluate_clustering.py specifically for a
# BINge paper analysis. It is intended to facilitate a clustering evaluation
# where the true clusters are derived from the human-mouse-rat-zebrafish
# homolog groups.

import os, argparse, sys
import networkx as nx
from sklearn.metrics.cluster import adjusted_rand_score, rand_score, \
    normalized_mutual_info_score, adjusted_mutual_info_score

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_ClustIO

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
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate input file format
    if args.clusterer == "binge":
        isBinge = validate_cluster_file(args.clusterFile)
        if not isBinge:
            raise ValueError("The input file was not validated as a BINge file!")
    elif args.clusterer == "cdhit":
        isBinge = validate_cluster_file(args.clusterFile)
        if isBinge:
            raise ValueError("The input file was not validated as a CD-HIT file!")
    elif args.clusterer == "corset" or args.clusterer == "mmseqs":
        isTSV = validate_cluster_tsv_file(args.clusterFile)
    else:
        raise NotImplementedError()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def validate_cluster_tsv_file(fileName):
    '''
    Validation function specifically for handling Corset-like or MMseqs2 cluster results.
    The file is expected to strictly conform to there being only two columns.
    
    Parameters:
        fileName -- a string indicating the location of the TSV file for validation.
    Returns:
        isValid -- a boolean that will be True; the alternative is quit() being called.
    '''
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l != "":
                sl = l.split("\t")
                if len(sl) != 2:
                    errorMsg = (f"The input file '{fileName}' does not appear to be a TSV (e.g., Corset) " +
                        "cluster file.\nYou should check your inputs and try again.")
                    raise ValueError(errorMsg)
    return True

def parse_corset_clusters(fileName):
    '''
    After a corset TSV has been validated, this function will parse it to a dictionary
    structure. The left column of the TSV must be the sequence ID, with the right
    column being the cluster ID.
    
    Parameters:
        fileName -- a string indicating the location of the Corset TSV file for parsing.
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                           0: ['seqID1', 'seqID2'],
                           1: ['seqID3'],
                           ...
                       }
    '''
    clusterDict = {}
    clusterNum = -1
    lastCluster = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l != "":
                seqID, clustID = l.split("\t")
                if clustID != lastCluster:
                    clusterNum += 1
                lastCluster = clustID
                
                clusterDict.setdefault(clusterNum, [])
                clusterDict[clusterNum].append(seqID)
    return clusterDict

def parse_mmseqs_clusters(fileName):
    '''
    After a MMseqs2 TSV has been validated, this function will parse it to a dictionary
    structure. The left column of the TSV must be the cluster representative ID, with the
    right column being the member sequence ID.
    
    Parameters:
        fileName -- a string indicating the location of the Corset TSV file for parsing.
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                           0: ['seqID1', 'seqID2'],
                           1: ['seqID3'],
                           ...
                       }
    '''
    clusterDict = {}
    clusterNum = -1
    lastCluster = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l != "":
                clustID, seqID = l.split("\t")
                if clustID != lastCluster:
                    clusterNum += 1
                lastCluster = clustID
                
                clusterDict.setdefault(clusterNum, [])
                clusterDict[clusterNum].append(seqID)
    return clusterDict

def parse_orthologs(orthologsFile):
    '''
    Parses a gene_orthologs file from NCBI to a dictionary structure.
    
    Parameters:
        orthologsFile -- a string indicating the location of the gene_orthologs file.
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
        
        # Handle HOM_AllOrganism.rpt file
        if firstLine.startswith("DB Class Key"):
            found = set()
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                dbClassKey, geneID = sl[0], sl[4]
                
                if not dbClassKey in found:
                    numGeneClusters += 1
                    found.add(dbClassKey)
                
                trueDict[geneID] = numGeneClusters
        
        # Handle gene_orthologs file
        elif firstLine.startswith("#tax_id\tGeneID"):
            graph = nx.Graph()
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                geneID, otherGeneID = sl[0], sl[1]
                
                graph.add_node(geneID)
                graph.add_node(otherGeneID)
                graph.add_edge(geneID, otherGeneID)
            
            for connectedSeqIDs in nx.connected_components(graph):
                numGeneClusters += 1
                for seqID in connectedSeqIDs:
                    trueDict[seqID] = numGeneClusters
        
        # Raise error for unhandled file
        else:
            raise ValueError("The input file does not appear to be a gene_orthologs or HOM_AllOrganism.rpt file!")
    
    return trueDict, numGeneClusters

def parse_fasta_geneids(refseqFastaFile):
    '''
    Parses a RefSeq annotations FASTA from NCBI to a dictionary structure.
    
    Parameters:
        refseqFastaFile -- a string indicating the location of the FASTA file.
    Returns:
        idMappingDict -- a dictionary with structure like:
                    {
                        'seqID1': 'entrezID1'
                        'seqID2': 'entrezID2',
                        ...
                    }
    '''
    idMappingDict = {}
    with open(refseqFastaFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                seqID = line[1:].split(" ")[0]
                dbxref = line.split("[db_xref=")[1].split("]")[0]
                
                geneID = None
                for xref in dbxref.split(","):
                    if xref.startswith("GeneID:"):
                        geneID = xref.split(":")[1]
                        break
                assert geneID != None
                
                idMappingDict[seqID] = geneID
    return idMappingDict

## Main
def main():
    # User input
    usage = """%(prog)s is a modified evaluate_clustering.py script designed to receive not a
    reference annotation file, but a gene_orthologs file from https://ftp.ncbi.nih.gov/gene/DATA
    or the HOM_AllOrganism.rpt file from https://www.informatics.jax.org.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-g", dest="geneOrthologsFile",
                   required=True,
                   help="Input gene_orthologs or HOM_AllOrganism.rpt file")
    p.add_argument("-c", dest="clusterFile",
                   required=True,
                   help="Input the CD-HIT or BINge cluster file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input the RefSeq FASTA file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for text results")
    p.add_argument("-p", dest="clusterer",
                   required=True,
                   choices=["binge", "cdhit", "corset", "mmseqs"],
                   help="Specify which clusterer's results you are providing.")
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
    p.add_argument("--noMRNA", dest="noMRNA",
                   required=False,
                   action="store_true",
                   help="""Optionally, if your genome is bacterial or archaeal, then your GFF3
                   likely does not have mRNA features; it lists CDS directly under the gene
                   feature. Specify this flag to allow for that behaviour ONLY if you are
                   looking at one of these organisms.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the orthologs file
    trueDict, numGeneClusters = parse_orthologs(args.geneOrthologsFile)
    
    # Parse the CD-HIT / BINge cluster file, changing cluster IDs to not overlap
    if args.clusterer == "binge":
        testDict = parse_binge_clusters(args.clusterFile) # we will 'test' this against our ground truth
    elif args.clusterer == "cdhit":
        testDict = ZS_ClustIO.CDHIT.parse_clstr_file(args.clusterFile) 
    elif args.clusterer == "corset":
        testDict = parse_corset_clusters(args.clusterFile)
    elif args.clusterer == "mmseqs":
        testDict = parse_mmseqs_clusters(args.clusterFile)
    else:
        raise NotImplementedError()
    
    # Parse the RefSeq file
    idMappingDict = parse_fasta_geneids(args.fastaFile)
    
    # Flip the dict around and +numGeneClusters to prevent cluster number overlap
    testDict = {
        seqID : clustNum+numGeneClusters
        for clustNum, idList in testDict.items()
        for seqID in idList
    }
    priorSize = len(testDict)
    
    # Modify RefSeq IDs to be entrez gene IDs
    testDict = {
        idMappingDict[seqID] : clustNum
        for seqID, clustNum in testDict.items()
        if seqID in idMappingDict # RefSeq CDS models have absences for some reason
    }
    postSize = len(testDict)
    print(f"Dropped {priorSize-postSize} sequences when getting entrez gene IDs")
    print(f"Final number of sequences to test: {postSize}")
    
    # Drop any sequences in testDict that aren't in our trueDict
    "Must have the exact same sequences for comparison"
    for seqID in list(testDict.keys()):
        if not seqID in trueDict:
            print(f"Dropped {seqID} from test labels as it has no match in true labels")
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
                    print(f"Dropped {seqID} from TRUE labels as it has no match in test labels")
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
    else:
        if args.noMRNA == True:
            print("I see {0} gene features in the GFF3".format(len(gff3.types["gene"])))
        else:
            print("I see {0} mRNA features in the GFF3".format(len(gff3.types["mRNA"])))
        print("From that, I've found {0} testable mRNA features".format(len(trueDict)))
        print(f"And I see {len(testDict)} sequences in the cluster file (after pruning to match the GFF3)")
        print("I will continue the evaluation, but if these numbers look odd then something may have gone wrong.")
    
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
        fileOut.write(f"Orthologs file\t{args.geneOrthologsFile}\n")
        fileOut.write(f"FASTA file\t{args.fastaFile}\n")
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
