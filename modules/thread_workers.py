import os, sys, time
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .bins import BinCollection, Bin, BinSplitter
from .gff3_handling import iterate_through_gff3
from threading import Thread

####

class GmapBinThread(Thread):
    '''
    This provides a modified Thread which allows for the output of bin_by_gmap
    to be stored locally within the object. Multi-threading this process is viable
    since each GMAP file should correspond to an entirely separate genome and hence
    there should be no overlap.
    '''
    def __init__(self, gmapFile, binCollection):
        Thread.__init__(self)
        
        self.gmapFile = gmapFile
        self.binCollection = binCollection
        self.novelBinCollection = None
        self.multiOverlaps = None
        self.exception = None
    
    def bin_by_gmap(self):
        '''
        Self parameters:
            gmapFile -- a string pointing to the location of a GMAP GFF3 file
                        containing transcript alignments.
            binCollection -- a BinCollection object containing gene feature bins
                            generated from an official genome annotation.
        '''
        # Behavioural parameters (static for now, may change later)
        OKAY_COVERAGE = 96.5
        OKAY_IDENTITY = 95
        OKAY_INDEL_PROPORTION = 0.01
        ##
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        ##
        GENE_BIN_OVL_PCT = 0.5
        NOVEL_BIN_OVL_PCT = 0.5
        
        # Parse GMAP file and begin binning
        novelBinCollection = BinCollection() # holds onto novel bins we create during this process
        pathDict = {} # holds onto sequences we've already checked a path for
        multiOverlaps = []
        
        for feature in iterate_through_gff3(self.gmapFile):
            mrnaFeature = feature.mRNA[0]
            
            # Get alignment statistics
            coverage, identity = float(mrnaFeature.coverage), float(mrnaFeature.identity)
            exonList = Bin.format_exons_from_gff3_feature(mrnaFeature)
            exonLength = sum([
                end - start + 1
                for start, end in exonList
            ])
            indelProportion = int(mrnaFeature.indels) / exonLength
            
            # Skip processing if the alignment sucks
            isGoodAlignment = coverage >= OKAY_COVERAGE \
                              and identity >= OKAY_IDENTITY \
                              and indelProportion <= OKAY_INDEL_PROPORTION
            if not isGoodAlignment:
                continue
            
            # See if this path should be skipped because another better one was processed
            baseID = feature.ID.rsplit(".", maxsplit=1)[0]
            if baseID in pathDict:
                prevCov, prevIdent = pathDict[baseID]
                
                # Skip if the first path was the best
                if (coverage + ALLOWED_COV_DIFF) < prevCov \
                    or (identity + ALLOWED_IDENT_DIFF) < prevIdent:
                    continue
            
            # Now that we've checked if this path is good, we'll store it
            """This means on the next loop if there's another path for this gene, 
            but it's worse than this current one, we'll just skip it"""
            pathDict.setdefault(baseID, [coverage, identity])
            
            # See if this overlaps an existing bin
            binOverlap = self.binCollection.find(feature.contig, feature.start, feature.end)
            
            # If it does not overlap an existing bin ...
            if len(binOverlap) == 0:
                # ... either add it to an existing novel bin or create a new novel bin
                _novel_binner(mrnaFeature, novelBinCollection, NOVEL_BIN_OVL_PCT) # pass mRNA
            
            # Exclude any transcripts that overlap multiple genes
            elif len(binOverlap) > 1:
                "We expect these occurrences to be chimeric transcripts which should be filtered"
                multiOverlaps.append(feature)
                continue
            
            # Compare to the existing bin
            else:
                geneBin = binOverlap[0]
                featureOvlPct, binOvlPct = _calculate_overlap_percentages(feature, geneBin)
                shouldJoin = featureOvlPct >= GENE_BIN_OVL_PCT or binOvlPct >= GENE_BIN_OVL_PCT
                
                # If this should join the gene bin, do so now
                if shouldJoin:
                    geneBin.add(baseID, exonList)
                
                # Otherwise, add it to an existing novel bin or create a new novel bin
                else:
                    _novel_binner(mrnaFeature, novelBinCollection, NOVEL_BIN_OVL_PCT)
        
        return novelBinCollection, multiOverlaps
    
    def run(self):
        try:
            self.novelBinCollection, self.multiOverlaps = self.bin_by_gmap()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception

def _novel_binner(mrnaFeature, novelBinCollection, NOVEL_BIN_OVL_PCT):
    '''
    Helper function pulled out to prevent code repetition. It will
    take a feature that has been determined to be "novel" and handle
    whether it should be binned with another novel bin or become
    a new bin itself.
    '''
    baseID = mrnaFeature.ID.rsplit(".", maxsplit=1)[0]
    
    thisBin = Bin(mrnaFeature.contig, mrnaFeature.start, mrnaFeature.end) # create a novel bin
    thisBin.add(baseID, Bin.format_exons_from_gff3_feature(mrnaFeature)) # we don't want the .path# / .mrna# suffix attached to this bin
    
    # ... see if this overlaps a novel bin that it should join
    novelBinOverlap = novelBinCollection.find(mrnaFeature.contig, mrnaFeature.start, mrnaFeature.end)
    binsToJoin = []
    for novelBin in novelBinOverlap:
        featureOvlPct, binOvlPct = _calculate_overlap_percentages(mrnaFeature, novelBin)
        shouldJoin = featureOvlPct >= NOVEL_BIN_OVL_PCT or binOvlPct >= NOVEL_BIN_OVL_PCT
        if shouldJoin:
            binsToJoin.append(novelBin)
    
    # Merge any bins that should join based on this novel sequence
    for novelBin in binsToJoin:
        thisBin.merge(novelBin)
        novelBinCollection.delete(novelBin)
    
    # Store this bin in the collection
    novelBinCollection.add(thisBin)

def _calculate_overlap_percentages(feature, bin):
    '''
    Feature should come from the GMAP parsing.
    Bin should come from a BinCollection.
    '''
    overlapLength = min(feature.end, bin.end) - max(feature.start, bin.start)
    featureOvlPct = overlapLength / (feature.end - feature.start)
    binOvlPct = overlapLength / (bin.end - bin.start)
    
    return featureOvlPct, binOvlPct

####

class BinSplitWorkerThread(Thread):
    '''
    This provides a modified Thread which is only being done for code tidiness and OOP
    reasons. Encapsulates the code needed to split a bin using the BinSplitter Class.
    
    Parameters:
        binQueue -- a queue.Queue object containing Bin objects for splitting.
        outputQueue -- a queue.Queue object that will receive the processed outputs
                       of this worker thread.
    '''
    def __init__(self, binQueue, outputQueue):
        Thread.__init__(self)
        
        self.binQueue = binQueue
        self.outputQueue = outputQueue
        self.exception = None
    
    def worker_function(self):
        '''
        Receives Bin objects via self.binQueue and processes them, adding
        the bin and the BinSplitter results to self.outputQueue.
        '''
        while True:
            # Continue condition
            if self.binQueue.empty():
                time.sleep(0.5)
                continue
            
            # Grabbing condition
            bin = self.binQueue.get()
            
            # Exit condition
            if bin == None:
                self.binQueue.task_done()
                break
            
            # Perform work
            splitter = BinSplitter(bin)
            splitter.cluster()
            clusterDict = splitter.resultClusters
            
            # Store result in output queue
            self.outputQueue.put([bin, clusterDict])
            
            # Mark work completion
            self.binQueue.task_done()
    
    def run(self):
        try:
            self.worker_function()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception

class CollectionWorkerThread(Thread):
    '''
    This provides a modified Thread which allows for multiple BinSplitter objects to
    run in parallel, feeding their output to an instance of this Class which combines
    their results into a new BinCollection.
    
    Parameters:
        outputQueue -- a queue.Queue object that receives outputs from worker threads.
    '''
    def __init__(self, outputQueue):
        Thread.__init__(self)
        
        self.outputQueue = outputQueue
        self.binCollection = None
        self.exception = None
    
    def worker_function(self):
        '''
        Receives Bin objects and their BinSplitter dictionary results via self.outputQueue
        and stores them in a new BinCollection.
        '''
        self.binCollection = BinCollection()
        
        # Process results as they become available
        while True:
            # Continue condition
            if self.outputQueue.empty():
                time.sleep(0.5)
                continue
            
            # Grabbing condition
            bin, clusterDict = self.outputQueue.get()
            
            # Exit condition
            if clusterDict == None:
                self.outputQueue.task_done()
                break
            
            # Perform work
            for clusterIDs in clusterDict.values():
                # Derive the position of this cluster
                contig = bin.contig
                start = min([ start for seqID in clusterIDs for start, end in bin.exons[seqID] ])
                end = max([ end for seqID in clusterIDs for start, end in bin.exons[seqID] ])
                exonsDict = { seqID : [] for seqID in clusterIDs } # exons are not needed anymore
                
                # Create a bin for this split cluster
                newBin = Bin(contig, start, end)
                newBin.union(clusterIDs, exonsDict)
                
                # Store it in the BinCollection
                self.binCollection.add(newBin)
            
            # Mark work completion
            self.outputQueue.task_done()
    
    def run(self):
        try:
            self.worker_function()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception
