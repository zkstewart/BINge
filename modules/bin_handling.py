import os, queue

from .gff3_handling import GFF3
from .thread_workers import GmapBinThread, BinSplitWorkerThread, CollectionWorkerThread
from .bins import Bin, BinCollection

def generate_bin_collections(workingDirectory):
    '''
    Receives a list of genome FASTA files and generates bin collections each of these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory with genome GFF3
                            and/or FASTA files in the subdirectory 'genomes'.
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
    
    # Generate and seed bin collections
    collectionList = []
    
    for gff3File, genomeFile, suffixNum in filePairs:
        binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if gff3File != None:
            gff3Obj = GFF3(gff3File, strict_parse=False)
            for geneFeature in gff3Obj.types["gene"]:
                
                # Create a bin for this feature
                featureBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
                try:
                    featureBin.add(geneFeature.ID, Bin.format_exons_from_gff3_feature(geneFeature))
                except:
                    "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                    continue
                
                # See if this overlaps an existing bin
                binOverlap = binCollection.find(geneFeature.contig, geneFeature.start, geneFeature.end)
                
                # If not, add the new bin
                if len(binOverlap) == 0:
                    binCollection.add(featureBin)
                
                # Otherwise...
                else:
                    # ... merge the bins together
                    for overlappingBin in binOverlap:
                        featureBin.merge(overlappingBin)
                    
                    # ... delete the overlapping bins
                    for overlappingBin in binOverlap:
                        binCollection.delete(overlappingBin)
                    
                    # ... and add the new bin
                    binCollection.add(featureBin)
        
        # Store the bin collection in our list for multi-threading later
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
        gmapIdentity -- a list of floats indicating what identity value a GMAP alignment
                        must have for it to be considered for binning.
    Modifies:
        collectionList -- the BinCollection objects in this list will have alignment IDs
                          added to the contained Bin objects.
    Returns:
        novelBinCollection -- a BinCollection containing Bins which did not overlap existing
                              Bins part of the input collectionList.
    '''
    novelBinCollection = BinCollection()
    multiOverlaps = []
    for i in range(0, len(gmapFiles), threads): # only process n (threads) files at a time
        processing = []
        for x in range(threads): # begin processing n files
            if i+x < len(gmapFiles): # parent loop may excess if n > the number of GMAP files
                gmapFile = gmapFiles[i+x]
                binCollection = collectionList[i+x]
                thisGmapIdentity = gmapIdentity[i+x]
                
                gmapWorkerThread = GmapBinThread(gmapFile, binCollection, thisGmapIdentity)
                processing.append(gmapWorkerThread)
                gmapWorkerThread.start()
        
        # Gather results
        for gmapWorkerThread in processing:
            # Wait for thread to end
            gmapWorkerThread.join()
            
            # Grab the new outputs from this thread
            "Each thread modifies a BinCollection part of collectionList directly"
            threadNovelBinCollection = gmapWorkerThread.novelBinCollection
            multiOverlaps.append(gmapWorkerThread.multiOverlaps)
            
            # Merge them via flattening
            novelBinCollection.flatten(threadNovelBinCollection)
    return novelBinCollection, multiOverlaps

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

def multithread_bin_splitter(binCollection, threads):
    '''
    Using a queue system, worker threads will grab bins from the binCollection parameter
    process with the BinSplitter class. They will put their results into the output
    worker which will help to format a new BinCollection object as output.
    
    We want to split bins because each bin is only guaranteed to contain
    sequences that generally align well over a predicted genomic region. There's no
    guarantee that they are isoforms of the same gene, or even share any exonic content.
    This splitting step is very lax and designed specifically for separating out nested
    genes within the introns of what might be considered the "main gene" the
    bin represents.
    
    Parameters:
        binCollection -- a BinCollection containing Bins which have aligned against
                         a reference genome, and hence have informative .exons values.
        threads -- an integer indicating how many threads to run.
        clusterType -- a string to put in an output column indicating what source of
                       evidence led to this cluster's formation; usually, this will
                       have the value of "binned" to indicate that it was binned via
                       GMAP alignment, or "unbinned" to indicate that it has been
                       binned solely via CD-HIT clustering.
    Returns:
        newBinCollection -- a new BinCollection object to replace the original one.
    '''
    # Set up queueing system for multi-threading
    workerQueue = queue.Queue(maxsize=int(threads * 10)) # allow queue to scale with threads
    outputQueue = queue.Queue(maxsize=int(threads * 20))
    
    # Start up threads for clustering of bins
    workers = []
    for _ in range(threads):
        worker = BinSplitWorkerThread(workerQueue, outputQueue)
        worker.daemon = True
        worker.start()
        workers.append(worker)
    
    outputWorker = CollectionWorkerThread(outputQueue)
    outputWorker.daemon = True
    outputWorker.start()
    
    # Put bins in queue for worker threads
    for interval in binCollection:
        bin = interval.data
        workerQueue.put(bin)
    
    # Close up shop on the threading structures
    for i in range(threads):
        workerQueue.put(None) # this is a marker for the worker threads to stop
    for worker in workers:
        worker.join() # need to call .join() on the workers themselves, not the queue!
    
    outputQueue.put([None, None]) # marker for the output thead to stop; needs 2 Nones!
    outputWorker.join()
    
    return outputWorker.binCollection
