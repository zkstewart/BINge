import os, pickle
from .salmon import EquivalenceClassCollection, QuantCollection, DGEQuantCollection

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

def parse_dge_quants(quantFiles, sampleNames):
    '''
    Parses in one or more quant files from Salmon, producing a DGEQuantCollection
    object enabling the generation of output files needed for DESeq2 DGE analysis.
    
    Parameters:
        quantFiles -- a list containing strings pointing to the location
                      of Salmon quant files (quant.sf files).
        sampleNames -- an equal length list indicating the sample names for each
                       salmon quant file.
    '''
    assert len(quantFiles) == len(sampleNames), \
        ("parse_dge_quants cannot parse quant files since the number " +
         f"of files ({len(quantFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    dgeQuantCollection = DGEQuantCollection()
    for i in range(len(quantFiles)):
        quantFile = quantFiles[i]
        sample = sampleNames[i]
        
        dgeQuantCollection.parse_quant_file(quantFile, sample)
    return dgeQuantCollection

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

def load_sequence_length_index(indexFile):
    '''
    Load in the pickled result of generate_sequence_length_index().
    
    Parameters:
        indexFile -- a string indicating the location of index generated from a FASTA file.
    '''
    if os.path.exists(indexFile):
        with open(indexFile, "rb") as fileIn:
            seqLenDict = pickle.load(fileIn)
        return seqLenDict
    else:
        raise FileNotFoundError(("load_sequence_length_index() failed because " + 
                                 f"'{indexFile}' doesn't exist!"))
