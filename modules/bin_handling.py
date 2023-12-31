import os
from multiprocessing import Queue

from .thread_workers import GmapBinProcess, FragmentFixProcess, \
    CollectionSeedProcess

TIMEOUT_LIMIT = 60 # seconds

def generate_bin_collections(workingDirectory, threads):
    '''
    Receives a list of genome FASTA files and generates bin collections each of these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory with genome GFF3
                            and/or FASTA files in the subdirectory 'genomes'.
        threads -- an integer indicating how many threads to run.
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
            filePairs.append([None, os.path.join(genomesDir, file), suffixNum])
            
            # Add any corresponding GFF3 file if one exists
            gff3File = f"annotation{suffixNum}.gff3"
            if os.path.exists(os.path.join(genomesDir, gff3File)):
                filePairs[-1][0] = os.path.join(genomesDir, gff3File)
    
    assert len(filePairs) > 0, \
        f"generate_bin_collections failed because '{genomesDir}' contains no genomes somehow?"
    
    # Sort the list for consistency of ordering
    "If there are gaps in the suffixNum's, ordering is important to keep things paired up"
    filePairs.sort(key = lambda x: int(x[2]))
    suffixes = [ int(x[2]) for x in filePairs ]
    isConsecutive = suffixes == list(range(1, len(suffixes)+1))
    assert isConsecutive, \
        (f"generate_bin_collections failed because genome files in '{genomesDir}' are not " + 
        "consecutively numbered, which is important for later program logic!")
    
    # Start up threads
    finished = []
    for i in range(0, len(filePairs), threads): # only process n (threads) at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(filePairs): # parent loop may excess if n > the number of GMAP files
                gff3File, _, _ = filePairs[i+x]
                
                seedWorkerThread = CollectionSeedProcess(gff3File)
                
                processing.append(seedWorkerThread)
                seedWorkerThread.start()
        
        # Wait on processes to end
        for seedWorkerThread in processing:
            seedWorkerThread.join()
            finished.append(seedWorkerThread)
    
    # Gather results
    collectionList = []
    for seedWorkerThread in finished:
        seedWorkerThread.check_errors()
        binCollection = seedWorkerThread.get_result(TIMEOUT_LIMIT)
        collectionList.append(binCollection)
    
    return collectionList

def populate_bin_collections(collectionList, gmapFiles, threads, gmapIdentity):
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
    Modifies:
        collectionList -- the BinCollection objects in this list will have alignment IDs
                          added to the contained Bin objects.
    Returns:
        novelBinCollection -- a BinCollection containing Bins which did not overlap existing
                              Bins part of the input collectionList.
    '''
    # Establish lists for feeding data into threads
    threadData = []
    for suffixNum in range(1, len(collectionList) + 1):
        # Figure out which genome we're looking at
        genomePrefix = f"genome{suffixNum}_"
        genomeIndex = int(suffixNum) - 1
        
        # Get all GMAP files associated with this genome
        thisGmapFiles = [ gmFile for gmFile in gmapFiles if genomePrefix in gmFile]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[genomeIndex]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection])
    
    # Start up threads
    finished = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection = threadData[i+x]
                
                populateWorkerThread = GmapBinProcess(thisGmapFiles, thisBinCollection,
                                                      gmapIdentity)
                
                processing.append(populateWorkerThread)
                populateWorkerThread.start()
        
        # Wait on processes to end
        for populateWorkerThread in processing:
            populateWorkerThread.join()
            finished.append(populateWorkerThread)
    
    # Gather results
    resultBinCollections, resultMultiOverlaps = [], []
    for populateWorkerThread in finished:
        populateWorkerThread.check_errors()
        binCollection, multiOverlap = populateWorkerThread.get_result(TIMEOUT_LIMIT)
        resultBinCollections.append(binCollection)
        resultMultiOverlaps.append(multiOverlap)
    
    return resultBinCollections, resultMultiOverlaps

def fix_collection_fragments(collectionList, multiOverlaps, threads):
    '''
    Connects bins which are putatively fragmented and should be merged
    into a single bin. The input collectionList and multiOverlaps objects
    should be ordered equivalently, such that each value in multiOverlaps
    is used for fragment fixing of its respectively indexed collectionList.
    
    Parameters:
        collectionList -- a list containing one or more BinCollection
                          objects.
        multiOverlaps -- a list containing GFF3 Feature objects indicating
                         GMAP alignments that overlapped more than one bin.
        threads -- an integer indicating how many threads to run.
    Returns:

    '''
    # Validate inputs
    assert len(collectionList) == len(multiOverlaps), \
        "fix_collection_fragments expects the two input lists to be equally sized!"
    assert isinstance(threads, int), \
        "fix_collection_fragments expects threads argument to be an integer!"
    assert threads >= 1, \
        "fix_collection_fragments expects threads argument to be >= 1"
    
    # Start up threads
    finished = []
    for i in range(0, len(collectionList), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(collectionList): # parent loop may excess if n > the number of GMAP files
                thisBinCollection, thisMultiOverlap = collectionList[i+x], multiOverlaps[i+x]
                
                fragmentWorkerThread = FragmentFixProcess(thisBinCollection,
                                                          thisMultiOverlap)
                
                processing.append(fragmentWorkerThread)
                fragmentWorkerThread.start()
        
        # Wait on processes to end
        for fragmentWorkerThread in processing:
            fragmentWorkerThread.join()
            finished.append(fragmentWorkerThread)
    
    # Gather results
    resultCollectionList = []
    for fragmentWorkerThread in finished:
        fragmentWorkerThread.check_errors()
        binCollection = fragmentWorkerThread.get_result(TIMEOUT_LIMIT)
        resultCollectionList.append(binCollection)
    
    return resultCollectionList

def iterative_bin_self_linking(binCollection, convergenceIters):
    '''
    Links the bins in a BinCollection to themselves iteratively until
    it converges or convergenceIters is reached.
    
    Parameters:
        binCollection -- a BinCollection object
        convergenceIters -- an integer providing a maximum limit for
                            the amount of iterations that may occur.
    '''
    prevCollectionCount = len(binCollection.bins)
    for _ in range(convergenceIters):
        binCollection = binCollection.link_bins()
        if len(binCollection.bins) == prevCollectionCount:
            break
        prevCollectionCount = len(binCollection.bins)
    return binCollection
