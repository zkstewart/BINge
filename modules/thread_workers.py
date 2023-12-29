import os, time, sys
from multiprocessing import Process, Pipe

from .gff3_handling import GFF3, iterate_through_gff3
from .bins import BinCollection, Bin, BinSplitter
from .bin_handling import add_bin_to_collection

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_MapIO import GMAP_DB

class GmapIndexThread(Process):
    '''
    This provides a modified Process which allows for GMAP indexing to be run
    in parallel using the GMAP_DB class.
    
    Parameters:
        fasta -- a string indicating the location of a FASTA file for GMAP indexing.
        gmapDir -- a string indicating the location of the GMAP executable files.
    '''
    def __init__(self, fasta, gmapDir):
        Process.__init__(self)
        
        # Processing defaults
        self.gmapFile = fasta
        self.gmapDir = gmapDir
        self.db = GMAP_DB(fasta, gmapDir)
        self.exception = None
    
    def generate_gmap_index(self):
        '''
        Self parameters:
            self.db -- a GMAP_DB object.
        '''
        if not self.db.index_exists():
            self.db.index()
    
    def run(self):
        try:
            self.generate_gmap_index()
        except BaseException as e:
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception

####

class GmapBinThread(Process):
    '''
    This provides a modified Process which allows for the output of bin_by_gmap
    to be stored locally within the object. Multi-processing this process is viable
    since each GMAP file should correspond to an entirely separate genome and hence
    there should be no overlap.
    
    Parameters:
        gmapFile -- a string indicating the location of a GMAP GFF3 for parsing.
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        connection -- a multiprocessing.Pipe() object to send the return value through.
        minIdentity -- a float fraction indicating what identity value is minimally required
                       for us to use a GMAP alignment.
    '''
    def __init__(self, gmapFiles, binCollection, multiOverlaps, connection, minIdentity=0.95):
        Process.__init__(self)
        
        # Behavioural parameters
        assert 0.0 < minIdentity <= 1.0, \
            "GmapBinThread needs to receive minIdentity as a float fraction from >0 to <=1"
        self.minIdentity = minIdentity * 100
        
        # Processing defaults
        self.gmapFiles = gmapFiles
        self.binCollection = binCollection
        self.multiOverlaps = multiOverlaps
        self.connection = connection
        
        # Exception handling
        self._exceptionReceiver, self._exceptionSender = Pipe()
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
        "These statics are for filtering GMAP alignments that are poor quality"
        OKAY_COVERAGE = 96.5
        OKAY_INDEL_PROPORTION = 0.01
        ##
        "These statics prevent a gene mapping to multiple locations when it has a clear best"
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        
        # Iterate through GMAP files
        for gmapFile in self.gmapFiles:
            pathDict = {} # holds onto sequences we've already checked a path for        
            for feature in iterate_through_gff3(gmapFile):
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
                                and identity >= self.minIdentity \
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
                
                # Now that we've checked if this path is good, we'll remember it
                """This means on the next loop if there's another path for this gene, 
                but it's worse than this current one, we'll just skip it"""
                pathDict.setdefault(baseID, [coverage, identity])
                
                # Create a new bin for this features
                newBin = _create_novel_bin(mrnaFeature, exonList)
                
                # See if this overlaps an existing bin
                coordinateOverlaps = self.binCollection.find(feature.contig, feature.start, feature.end)
                binOverlap = []
                for bin in coordinateOverlaps:
                    shouldMerge = newBin.is_overlapping(bin)
                    if shouldMerge:
                        binOverlap.append(bin)
                
                # If it does not overlap an existing bin ...
                if len(binOverlap) == 0:
                    # ... add the novel bin into the collection
                    self.binCollection.add(newBin)
                
                # Exclude any transcripts that overlap multiple bins
                elif len(binOverlap) > 1:
                    """We expect these occurrences to be chimeric transcripts which should
                    be filtered. Otherwise, fragmented transcripts being loaded in may
                    result in a single loci having multiple bins over it if the full length
                    one doesn't get processed first. We'll need to fix those later by
                    assessing these multi-overlaps."""
                    self.multiOverlaps.append(feature)
                
                # If it overlaps exactly 1 bin ...
                else:
                    # ... merge the bins together
                    add_bin_to_collection(self.binCollection, binOverlap, newBin)
    
    def run(self):
        try:
            self.bin_by_gmap()
            self.connection.send([self.binCollection, self.multiOverlaps])
        except BaseException as e:
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception

def _create_novel_bin(mrnaFeature, exonList):
    '''
    Helper function to create a Bin from a mRNA feature parsed out of a GFF3 file.
    exonList is fed in directly since it's already created in the parent function
    and this should offer a slight optimisation benefit.
    '''
    baseID = mrnaFeature.ID.rsplit(".", maxsplit=1)[0]
    thisBin = Bin(mrnaFeature.contig, mrnaFeature.start, mrnaFeature.end)
    thisBin.add(baseID, exonList) # we don't want the .path# / .mrna# suffix attached to this bin
    
    return thisBin

####

class BinSplitWorkerThread(Process):
    '''
    This provides a modified Process which is only being done for code tidiness and OOP
    reasons. Encapsulates the code needed to split a bin using the BinSplitter Class.
    
    Parameters:
        binQueue -- a multiprocess.Queue object containing Bin objects for splitting.
        outputQueue -- a multiprocess.Queue object that will receive the processed outputs
                       of this worker process.
    '''
    def __init__(self, binQueue, outputQueue):
        Process.__init__(self)
        
        self.binQueue = binQueue
        self.outputQueue = outputQueue
        
        # Exception handling
        self._exceptionReceiver, self._exceptionSender = Pipe()
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
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception

class CollectionWorkerThread(Process):
    '''
    This provides a modified Process which allows for multiple BinSplitter objects to
    run in parallel, feeding their output to an instance of this Class which combines
    their results into a new BinCollection.
    
    Parameters:
        outputQueue -- a multiprocess.Queue object that receives outputs from worker processes.
        connection -- a multiprocessing.Pipe() object to send the return value through.
    '''
    def __init__(self, outputQueue, connection):
        Process.__init__(self)
        
        self.outputQueue = outputQueue
        self.binCollection = None
        self.connection = connection
        
        # Exception handling
        self._exceptionReceiver, self._exceptionSender = Pipe()
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
            self.connection.send(self.binCollection)
        except BaseException as e:
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception

####

class FragmentFixThread(Process):
    '''
    This provides a modified Process which allows for fragment fixing to be run in
    parallel across multiple separate BinCollection() objects.
    
    Parameters:
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        multiOverlaps -- a list containing GFF3 Feature objects indicating
                         GMAP alignments that overlapped more than one bin.
        connection -- a multiprocessing.Pipe() object to send the return value through.
    '''
    def __init__(self, binCollection, multiOverlaps, connection):
        Process.__init__(self)
        
        self.binCollection = binCollection
        self.multiOverlaps = multiOverlaps
        self.connection = connection
        
        # Exception handling
        self._exceptionReceiver, self._exceptionSender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            self.binCollection.fix_fragments(self.multiOverlaps)
            self.connection.send(self.binCollection)
        except BaseException as e:
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception

####

class CollectionSeedThread(Process):
    '''
    This provides a modified Process which allows for generate_bin_collections() to
    be run in parallel for multiple genomes.
    
    Parameters:
        gff3File -- a string pointing to a GFF3 file giving annotations for a genome,
                    or None if no pre-seeding is to occur.
        connection -- a multiprocessing.Pipe() object to send the return value through.
    '''
    def __init__(self, gff3File, connection):
        Process.__init__(self)
        
        self.gff3File = gff3File
        self.binCollection = None
        self.connection = connection
        
        # Exception handling
        self._exceptionReceiver, self._exceptionSender = Pipe()
        self.exception = None
    
    def bin_seeder(self):
        self.binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if self.gff3File != None:
            gff3Obj = GFF3(self.gff3File, strict_parse=False)
            for geneFeature in gff3Obj.types["gene"]:
                
                # Create a bin for this feature
                featureBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
                try:
                    for mrnaFeature in geneFeature.mRNA:
                        featureBin.add(mrnaFeature.ID, Bin.format_exons_from_gff3_feature(mrnaFeature))
                except:
                    "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                    continue
                
                # See if this overlaps an existing bin
                binOverlap = self.binCollection.find(geneFeature.contig, geneFeature.start, geneFeature.end)
                
                # Handle the storing of this bin in its collection
                add_bin_to_collection(self.binCollection, binOverlap, featureBin)
    
    def run(self):
        try:
            self.bin_seeder()
            self.connection.send(self.binCollection)
        except BaseException as e:
            self._exceptionSender.send(e)
    
    def join(self):
        Process.join(self)
        if self._exceptionReceiver.poll():
            self.exception = self._exceptionReceiver.recv()
        if self.exception:
            raise self.exception
