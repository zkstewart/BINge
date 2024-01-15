import os

from .thread_workers import CollectionSeedProcess, GmapBinProcess, GraphPruneProcess

def generate_bin_collections(workingDirectory, threads, isMicrobial):
    '''
    Receives a list of genome FASTA files and generates bin collections each of these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory with genome GFF3
                            and/or FASTA files in the subdirectory 'genomes'.
        threads -- an integer indicating how many threads to run.
        isMicrobial -- a boolean indicating whether the genomes are microbial or not which,
                       in turn, determines whether we will parse mRNA features (False) or
                       gene features (True).
    Returns:
        collectionList -- a list containing BinCollections in numerical order of the genome
                          files.
    '''
    # Locate subdirectory containing files
    genomesDir = os.path.join(workingDirectory, "genomes")
    assert os.path.isdir(genomesDir), \
        f"generate_bin_collections failed because '{genomesDir}' isn't a directory somehow?"
    
    # Locate all genome/GFF3 pairings
    filePairs = []
    for file in os.listdir(genomesDir):
        if file.endswith(".fasta"):
            assert file.startswith("genome"), \
                f"FASTA file in '{genomesDir}' has a different name than expected?"
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".fasta")[0]
            suffixNum = filePrefix.split("genome")[1]
            assert suffixNum.isdigit(), \
                f"FASTA file in '{genomesDir}' does not have a number suffix?"
            
            # Add value to pairings list
            filePairs.append([None, suffixNum])
            
            # Add any corresponding GFF3 file if one exists
            gff3File = f"annotation{suffixNum}.gff3"
            if os.path.exists(os.path.join(genomesDir, gff3File)):
                filePairs[-1][0] = os.path.join(genomesDir, gff3File)
    
    assert len(filePairs) > 0, \
        f"generate_bin_collections failed because '{genomesDir}' contains no genomes somehow?"
    
    # Sort the list for consistency of ordering
    "If there are gaps in the suffixNum's, ordering is important to keep things paired up"
    filePairs.sort(key = lambda x: int(x[1]))
    suffixes = [ int(x[1]) for x in filePairs ]
    isConsecutive = suffixes == list(range(1, len(suffixes)+1))
    assert isConsecutive, \
        (f"generate_bin_collections failed because genome files in '{genomesDir}' are not " + 
        "consecutively numbered, which is important for later program logic!")
    
    # Start up threads
    collectionList = []
    for i in range(0, len(filePairs), threads): # only process n (threads) at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(filePairs): # parent loop may excess if n > the number of GMAP files
                gff3File, _ = filePairs[i+x]
                
                seedWorkerThread = CollectionSeedProcess(gff3File, isMicrobial)
                seedWorkerThread.start()
                processing.append(seedWorkerThread)
        
        # Wait on processes to end
        for seedWorkerThread in processing:
            binCollection = seedWorkerThread.get_result()
            seedWorkerThread.join()
            seedWorkerThread.check_errors()
            
            collectionList.append(binCollection)
    
    return collectionList

def populate_bin_collections(workingDirectory, collectionList, gmapFiles,
                             threads, gmapIdentity):
    '''
    Receives a list of BinCollection objects, alongside a list of GMAP GFF3 file
    locations, and uses multiple threads to parse the GMAP files and add them into
    an existing bin, or into a novel bin collection which is returned.
    
    Parameters:
        collectionList -- a list of BinCollection objects as resulting from
                          generate_bin_collections().
        gmapFiles -- a list of strings pointing to GMAP GFF3 files; they should
                     be ordered the same as the files which were used as input to
                     generate_bin_collections().
        threads -- an integer indicating how many threads to run; this code is
                   parallelised in terms of processing multiple GMAP files at a time,
                   if you have only 1 GMAP file then only 1 thread will be used.
        gmapIdentity -- a float indicating what identity value a GMAP alignment
                        must have for it to be considered for binning.
    Returns:
        novelBinCollection -- a BinCollection containing Bins which did not overlap existing
                              Bins part of the input collectionList.
    '''
    # Establish lists for feeding data into threads
    threadData = []
    for suffixNum in range(1, len(collectionList) + 1):
        # Figure out which genome we're looking at
        genomePrefix = f"genome{suffixNum}"
        genomeIndex = int(suffixNum) - 1
        
        # Get the location of the genome length index
        thisGenomeIndex = os.path.join(workingDirectory, "genomes",
                                       f"{genomePrefix}.fasta.lengths.pkl")
        
        # Get all GMAP files associated with this genome
        thisGmapFiles = [ gmFile for gmFile in gmapFiles if f"{genomePrefix}_" in gmFile]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[genomeIndex]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection, thisGenomeIndex])
    
    # Start up threads
    resultBinCollections = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection, thisGenomeIndex = threadData[i+x]
                
                populateWorkerThread = GmapBinProcess(thisGmapFiles, thisBinCollection,
                                                      thisGenomeIndex, gmapIdentity)
                populateWorkerThread.start()
                processing.append(populateWorkerThread)
        
        # Wait on processes to end
        for populateWorkerThread in processing:
            binCollection = populateWorkerThread.get_result()
            populateWorkerThread.join()
            populateWorkerThread.check_errors()
            
            resultBinCollections.append(binCollection)
        
    return resultBinCollections

def prune_graphs(binGraphList, threads):
    '''
    Receives a list of BinGraph objects and calls their .prune() method in parallel.
    
    Parameters:
        binGraphList -- a list containing BinGraph objects.
        threads -- an integer indicating how many threads to run.
    Returns:
        binGraphList -- a list containing the same BinGraph objects, just pruned.
    '''
    outputChimers = set()
    for i in range(0, len(binGraphList), threads): # only process n (threads) at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(binGraphList): # parent loop may excess if n > the number of GMAP files
                thisBinGraph = binGraphList[i+x]
                
                pruneWorkerThread = GraphPruneProcess(thisBinGraph)
                pruneWorkerThread.start()
                processing.append(pruneWorkerThread)
        
        # Wait on processes to end
        for pruneWorkerThread in processing:
            chimers = pruneWorkerThread.get_result()
            pruneWorkerThread.join()
            pruneWorkerThread.check_errors()
            
            outputChimers = outputChimers.union(chimers)
    
    return outputChimers
