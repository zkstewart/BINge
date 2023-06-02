import os, sys, time
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts import ZS_SeqIO, ZS_ClustIO

from .bins import BinCollection, Bin
from .gff3_handling import iterate_through_gff3
from threading import Thread

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
    
    def bin_by_gmap(self):
        '''
        Self parameters:
            gmapFile -- a string pointing to the location of a GMAP GFF3 file
                        containing transcript alignments.
            binCollection -- a BinCollection object containing gene feature bins
                            generated from an official genome annotation.
        '''
        # Behavioural parameters (static for now, may change later)
        OKAY_COVERAGE = 94
        OKAY_IDENTITY = 92
        ##
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        ##
        GENE_BIN_OVL_PCT = 0.5
        NOVEL_BIN_OVL_PCT = 0.5
        
        # Parse GMAP file and begin binning
        novelBinCollection = BinCollection() # holds onto novel bins we create during this process
        pathDict = {} # holds onto sequences we've already checked a path for
            
        for feature in iterate_through_gff3(self.gmapFile):
            mrnaFeature = feature.mRNA[0]
            
            # Get alignment statistics
            coverage, identity = float(mrnaFeature.coverage), float(mrnaFeature.identity)
            
            # Skip processing if the alignment sucks
            isGoodAlignment = coverage >= OKAY_COVERAGE and identity >= OKAY_IDENTITY
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
            """This means on the next loop if there's another path for this gene, but it's worse
            than this current one, we'll just skip it"""
            pathDict.setdefault(baseID, [coverage, identity])
            
            # See if this overlaps an existing bin
            binOverlap = self.binCollection.find(feature.contig, feature.start, feature.end)
            
            # If it does not overlap an existing bin ...
            if len(binOverlap) == 0:
                # ... either add it to an existing novel bin or create a new novel bin
                _novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT)
            
            # Exclude any transcripts that overlap multiple genes
            elif len(binOverlap) > 1:
                "We expect these occurrences to be chimeric transcripts which should be filtered"
                continue
            
            # Compare to the existing bin
            else:
                geneBin = binOverlap[0]
                featureOvlPct, binOvlPct = _calculate_overlap_percentages(feature, geneBin)
                shouldJoin = featureOvlPct >= GENE_BIN_OVL_PCT or binOvlPct >= GENE_BIN_OVL_PCT
                
                # If this should join the gene bin, do so now
                if shouldJoin:
                    geneBin.add(baseID)
                
                # Otherwise, add it to an existing novel bin or create a new novel bin
                else:
                    _novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT)
        
        return novelBinCollection
    
    def run(self):
        try:
            self.novelBinCollection = self.bin_by_gmap()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception

class OutputWorkerThread(Thread):
    '''
    This provides a modified Thread which allows for the output worker to return a value
    indicating how many clusters it wrote to file. This will save a little bit of time later
    when we want to append more cluster results to the same file.
    
    Parameters:
        outputQueue -- a queue.Queue object that receives outputs from worker threads.
        outputFileName -- a string indicating the location to write results to.
        clusterType -- a string to put in an output column indicating what source of
                       evidence led to this cluster's formation; usually, this will
                       have the value of "binned" to indicate that it was binned via
                       GMAP alignment, or "unbinned" to indicate that it has been
                       binned solely via CD-HIT clustering.
    '''
    def __init__(self, outputQueue, outputFileName, clusterType):
        Thread.__init__(self)
        
        self.outputQueue = outputQueue
        self.outputFileName = outputFileName
        self.clusterType = clusterType
        self.numClusters = None
    
    def worker_function(self):
        '''
        Receives Bin objects via self.outputQueue and writes them to file.
        Iterates self.numClusters to keep tally of how many clusters we've
        written. Because of the nature of this process, the value at the
        point of thread closure will be +1 higher than the actual number
        of written clusters.
        '''
        self.numClusters = 1
        with open(self.outputFileName, "w") as fileOut:
            # Write header
            fileOut.write("#BINge clustering information file\n")
            fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
            
            # Write content lines as they become available
            while True:
                # Continue condition
                if self.outputQueue.empty():
                    time.sleep(0.5)
                    continue
                
                # Grabbing condition
                clusterDict = self.outputQueue.get()
                
                # Exit condition
                if clusterDict == None:
                    self.outputQueue.task_done()
                    break
                
                # Perform work
                for clusterIDs in clusterDict.values():
                    for seqID in clusterIDs:
                        fileOut.write(f"{self.numClusters}\t{seqID}\t{self.clusterType}\n")
                    self.numClusters += 1
                
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

class WorkerThread(Thread):
    '''
    This provides a modified Thread which is only being done for code tidiness and OOP
    reasons. Encapsulates the code needed to cluster a Bin using CD-HIT.
    
    Parameters:
        binQueue -- a queue.Queue object containing Bin objects for clustering.
        outputQueue -- a queue.Queue object that will receive the processed outputs
                       of this worker thread.
        transcriptRecords -- a pyfaidx Fasta object containing all the sequences
                             we'll want to cluster from our bins.
        mem -- an integer indicating how many megabytes of memory to run CD-HIT with.
    '''
    def __init__(self, binQueue, outputQueue, transcriptRecords, mem):
        Thread.__init__(self)
        
        self.binQueue = binQueue
        self.outputQueue = outputQueue
        self.transcriptRecords = transcriptRecords
        self.mem = mem
    
    def worker_function(self):
        '''
        Receives Bin objects via self.binQueue and processes them, adding
        the results in self.outputQueue.
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
            clusterDict = self.cluster_bin(bin)
            
            # Store result in output queue
            self.outputQueue.put(clusterDict)
            
            # Mark work completion
            self.binQueue.task_done()
    
    def cluster_bin(self, bin):
        '''
        Performs the work of clustering the sequences in a bin with CD-HIT.
        
        Parameters:
            bin -- a Bin object which should be clustered with CD-HIT.
        
        Returns:
            clusterDict -- a dictionary result of clustering with structure like:
                        {
                            0: [seqid1, seqid2, ...],
                            1: [ ... ],
                            ...
                        }
        '''
        # Bypass clustering if this bin has only one sequence
        if len(bin.ids) == 1:
            clusterDict = {0: list(bin.ids)}
        else:
            # Create a FASTA object for the sequences in this bin
            FASTA_obj = ZS_SeqIO.FASTA(None)
            for id in bin.ids:
                FastASeq_obj = ZS_SeqIO.FastASeq(id, str(self.transcriptRecords[id]))
                FASTA_obj.insert(0, FastASeq_obj)
            
            # Cluster it with CD-HIT at lax settings
            clusterer = ZS_ClustIO.CDHIT(FASTA_obj, "nucleotide")
            clusterer.identity = 0.8
            clusterer.set_shorter_cov_pct(0.2)
            clusterer.set_longer_cov_pct(0.0)
            clusterer.set_local()
            clusterer.mem = self.mem
            clusterer.get_cdhit_results(returnFASTA=False, returnClusters=True) # throw away the temporary file return value
            
            clusterDict = clusterer.resultClusters
        
        return clusterDict
    
    def run(self):
        try:
            self.worker_function()
        except BaseException as e:
            self.exception = e
    
    def join(self):
        Thread.join(self)
        if self.exception:
            raise self.exception

def _novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT):
    '''
    Helper function pulled out to prevent code repetition. It will
    take a feature that has been determined to be "novel" and handle
    whether it should be binned with another novel bin or become
    a new bin itself.
    '''
    baseID = feature.ID.rsplit(".", maxsplit=1)[0]
    
    thisBin = Bin(feature.contig, feature.start, feature.end) # create a novel bin
    thisBin.add(baseID) # we don't want the .path# suffix attached to this bin
    
    # ... see if this overlaps a novel bin that it should join
    novelBinOverlap = novelBinCollection.find(feature.contig, feature.start, feature.end)
    binsToJoin = []
    for novelBin in novelBinOverlap:
        featureOvlPct, binOvlPct = _calculate_overlap_percentages(feature, novelBin)
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
    featureOvlPct = overlapLength / (feature.end - feature.start + 1)
    binOvlPct = overlapLength / (bin.end - bin.start + 1)
    
    return featureOvlPct, binOvlPct
