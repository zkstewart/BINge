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

def parse_gff3_ids(gff3File):
    '''
    Simple function to iterate through a GFF3 file and hold onto any parent and
    subfeature IDs. This process should encompass genes and mRNAs as well as any
    other subfeature e.g., lnc_RNA, whilst skipping over exons and CDS'.
    
    Parameters:
        gff3File -- a string indicating the location of a GFF3 file.
    Returns:
        annotIDs -- a set of strings containing the IDs of all parent and subfeatures
                    found in the GFF3 file.
    '''
    # Setup for slim parsing of attributes
    def _format_attributes(attributes):
        "Code borrowed from ZS_GFF3IO"
        SLIM_ATTRIBUTES = ["id", "parent"]
        
        splitAttributes = []
        for a in attributes.split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        attributesDict = {
            splitAttributes[i]: splitAttributes[i+1]
            for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)
            if splitAttributes[i].lower() in SLIM_ATTRIBUTES
        }
        return attributesDict
    
    parentIDs = set()
    subIDs = set()
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Parse content lines
            if not line.startswith("#") and len(sl) == 9:
                featureType, attributes = sl[2], sl[8]
                attributesDict = _format_attributes(attributes)
                if featureType == "gene":
                    parentIDs.add(attributesDict["ID"])
                elif "Parent" in attributesDict and attributesDict["Parent"] in parentIDs:
                    subIDs.add(attributesDict["ID"])
    annotIDs = parentIDs.union(subIDs)
    return annotIDs

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

def locate_read_files(readsDir, readsSuffix, isSingleEnd):
    '''
    Parameters:
        readsDir -- a list of strings indicating the directory where
                    reads files are located.
        readsSuffix -- a string indicating the suffix of the reads files.
        isSingleEnd -- a boolean indicating whether the reads are single-end (True)
                       or paired-end (False).
    Returns:
        forwardReads -- a list of strings indicating the location of forward reads.
        reverseReads -- a list of strings indicating the location of reverse reads
                        OR None if the reads are single-end.
        sampleNames -- a list of strings indicating the common prefix of the reads files
                       (if paired-end) or the file name sans suffix (if single-end).
    '''
    # Locate files from the directory
    forwardReads = []
    reverseReads = []
    sampleNames = []
    for dir in readsDir:
        for file in os.listdir(dir):
            if file.endswith(readsSuffix):
                if isSingleEnd:
                    forwardReads.append(os.path.join(dir, file))
                    sampleNames.append(file.rsplit(readsSuffix, maxsplit=1)[0])
                else:
                    if file.endswith(f"1{readsSuffix}"):
                        forwardReads.append(os.path.join(dir, file))
                    elif file.endswith(f"2{readsSuffix}"):
                        reverseReads.append(os.path.join(dir, file))
                    else:
                        raise ValueError(f"{file} ends with the expected suffix '{readsSuffix}' but is " +
                                        "not preceeded by a 1 or 2!")
    forwardReads.sort()
    reverseReads.sort()
    
    # Validate that paired files match
    if not isSingleEnd:
        if not len(forwardReads) == len(reverseReads):
            raise FileNotFoundError(f"Number of reads don't match for forward ({len(forwardReads)}) " + 
                                    f"and reverse ({len(reverseReads)}) files")
        
        for i in range(len(forwardReads)):
            forward = os.path.basename(forwardReads[i]) # strip directory which may falsely give a common prefix
            reverse = os.path.basename(reverseReads[i])
            
            prefix = os.path.commonprefix([forward, reverse]).rstrip("_ ")
            if prefix == "":
                raise ValueError("Forward and reverse read pairs don't have a common prefix after " +
                                 "sorting; make sure that files are named in a way that they can be " +
                                 "paired together!")
            if not forwardReads[i].startswith(f"{prefix}1") and reverseReads[i].startswith(f"{prefix}2"):
                raise ValueError(f"Forward read '{forwardReads[i]}' and reverse read '{forwardReads[i]}' " +
                                 "files don't start with a valid common prefix i.e., " + 
                                 f"'{prefix} should preceed a 1 or 2")
            sampleNames.append(prefix)
    
    # Return files
    return forwardReads, reverseReads if reverseReads != [] else None, sampleNames
