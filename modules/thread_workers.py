import os, sys, time
import networkx as nx
from multiprocessing import Process, Pipe, Queue

from .gff3_handling import GFF3, iterate_gmap_gff3
from .bins import BinCollection, Bin

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_MapIO import GMAP_DB

# Define helper functions
def add_bin_to_collection(binCollection, binOverlap, newBin):
    '''
    Assistant function to add a new bin on the basis of whatever bins it is
    overlapping. Modifies the binCollection input in place.
    
    Parameters:
        binCollection -- a BinCollection() object from which binOverlap was generated,
                         and the newBin is intended to be added to.
        binOverlap -- a list containing zero or more Bin() objects that the newBin
                      is overlapping.
        newBin -- a Bin() for addition to the binCollection input.
    '''
    # If the newBin does not overlap anything, add it in as-is
    if len(binOverlap) == 0:
        binCollection.add(newBin)
    
    # Otherwise...
    else:
        # ... merge any overlapping bins together
        for overlappingBin in binOverlap:
            newBin.merge(overlappingBin)
        
        # ... delete the overlapping bins
        for overlappingBin in binOverlap:
            binCollection.delete(overlappingBin)
        
        # ... and add the new bin in its place
        binCollection.add(newBin)

def _create_novel_bin(dataDict):
    thisBin = Bin(dataDict["contig"], dataDict["start"], dataDict["end"])
    thisBin.add(dataDict["Name"], dataDict["exons"])
    return thisBin

# Create base classes to inherit from
class BasicProcess(Process):
    def __init__(self, *args, **kwargs):
        Process.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self._exception_receiver, self._exception_sender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            self.task(*self.args, **self.kwargs)
        except Exception as e:
            self._exception_sender.send(e)
    
    def task(self, *args, **kwargs):
        # Override this method in a subclass
        raise NotImplementedError
    
    def check_errors(self):
        if self._exception_receiver.poll():
            self.exception = self._exception_receiver.recv()
        if self.exception:
            raise self.exception

class ReturningProcess(BasicProcess):
    def __init__(self, *args, **kwargs):
        BasicProcess.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.queue = Queue()
        self._exception_receiver, self._exception_sender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            result = self.task(*self.args, **self.kwargs)
            self.queue.put(result)
        except Exception as e:
            self._exception_sender.send(e)
            self.queue.put(None)
    
    def task(self, *args, **kwargs):
        # Override this method in a subclass
        raise NotImplementedError
    
    def get_result(self, timeout=None):
        return self.queue.get(timeout=timeout)

# Create inheriting classes
class GmapIndexProcess(BasicProcess):
    '''
    Handles GMAP indexing in a separate thread.
    
    Parameters:
        fasta -- a string indicating the location of a FASTA file for GMAP indexing.
        gmapDir -- a string indicating the location of the GMAP executable files.
    '''
    def task(self, fasta, gmapDir):
        db = GMAP_DB(fasta, gmapDir)
        if not db.index_exists():
            db.index()

class CollectionSeedProcess(ReturningProcess):
    '''
    Allows for generate_bin_collections() to be run in parallel for multiple genomes.
    
    Parameters:
        gff3File -- a string pointing to a GFF3 file giving annotations for a genome,
                    or None if no pre-seeding is to occur.
        queue -- a multiprocessing.Queue() object to send the return value through.
    '''
    def task(self, gff3File):
        binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if gff3File != None:
            gff3Obj = GFF3(gff3File, strict_parse=False)
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
                binOverlap = binCollection.find(geneFeature.contig, geneFeature.start, geneFeature.end)
                
                # Handle the storing of this bin in its collection
                add_bin_to_collection(binCollection, binOverlap, featureBin)
        
        return binCollection

class GmapBinProcess(ReturningProcess):
    '''
    Bins GMAP alignments into each genome's BinCollection object.
    
    Parameters:
        gmapFile -- a string indicating the location of a GMAP GFF3 for parsing.
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        minIdentity -- (optional) a float fraction indicating what identity value is
                       minimally required for us to use a GMAP alignment; default=0.95.
    '''
    def task(self, gmapFiles, binCollection, minIdentity=0.95):
        multiOverlaps = []
        
        # Behavioural parameters (static for now, may change later)
        "These statics are for filtering GMAP alignments that are poor quality"
        OKAY_COVERAGE = 96.5
        OKAY_INDEL_PROPORTION = 0.01
        ##
        "These statics prevent a gene mapping to multiple locations when it has a clear best"
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        
        # Iterate through GMAP files
        for gmapFile in gmapFiles:
            pathDict = {} # holds onto sequences we've already checked a path for        
            for dataDict in iterate_gmap_gff3(gmapFile):
                # Get alignment statistics
                coverage, identity = float(dataDict["coverage"]), float(dataDict["identity"])
                exonLength = sum([
                    end - start + 1
                    for start, end in dataDict["exons"]
                ])
                indelProportion = int(dataDict["indels"]) / exonLength
                
                # Skip processing if the alignment sucks
                isGoodAlignment = coverage >= OKAY_COVERAGE \
                                and identity >= minIdentity \
                                and indelProportion <= OKAY_INDEL_PROPORTION
                if not isGoodAlignment:
                    continue
                
                # See if this path should be skipped because another better one was processed
                if dataDict["Name"] in pathDict:
                    prevCov, prevIdent = pathDict[dataDict["Name"]]
                    
                    # Skip if the first path was the best
                    if (coverage + ALLOWED_COV_DIFF) < prevCov \
                        or (identity + ALLOWED_IDENT_DIFF) < prevIdent:
                        continue
                
                # Now that we've checked if this path is good, we'll remember it
                """This means on the next loop if there's another path for this gene, 
                but it's worse than this current one, we'll just skip it"""
                pathDict.setdefault(dataDict["Name"], [coverage, identity])
                
                # Create a new bin for this features
                newBin = _create_novel_bin(dataDict)
                
                # See if this overlaps an existing bin
                binOverlap = binCollection.find(dataDict["contig"], dataDict["start"],
                                                dataDict["end"])
                
                # Add a new bin, or merge any bins as appropriate
                add_bin_to_collection(binCollection, binOverlap, newBin)
        
        return binCollection, multiOverlaps

class QueuedBinSplitterProcess(ReturningProcess):
    '''
    Allows for bin splitting to be run in parallel across multiple separate
    BinCollection() objects. Utilises a queue to receive and send Bin() objects.
    
    Parameters:
        inputQueue -- a multiprocessing.Queue() object to receive Bin() objects from.
        shorterCovPct -- a float value indicating the percentage of the shortest
                         sequence's coverage that must be covered by the longer sequence;
                         default==0.40.
        longerCovPct -- a float value indicating the percentage of the longest
                        sequence's coverage that must be covered by the shorter sequence;
                        default==0.20.
    '''
    def task(self, inputQueue, shorterCovPct=0.40, longerCovPct=0.20):
        binCollection = BinCollection()
        
        while True:
            # Continue condition
            if inputQueue.empty():
                time.sleep(0.5)
                continue
            
            # Grabbing condition
            bin = inputQueue.get()
            
            # Exit condition
            if bin == None:
                break
            
            # Perform work
            binList = self._work_function(bin, shorterCovPct, longerCovPct)
            
            # Store results
            for bin in binList:
                binCollection.add(bin)
        
        return binCollection
    
    def _work_function(self, bin, shorterCovPct, longerCovPct):
        ids = list(bin.ids) # ensure we always use this in consistent ordering
        
        # Model links between sequences as a graph structure
        idsGraph = nx.Graph()
        idsGraph.add_nodes_from(ids)
        
        # Figure out which sequences have links through shared exons
        for i in range(len(ids)-1):
            for x in range(i+1, len(ids)):
                id1, id2 = ids[i], ids[x]
                exons1, exons2 = bin.exons[id1], bin.exons[id2]
                
                # Calculate the amount of overlap between their exons
                exonOverlap = sum([
                    abs(max(exons1[n][0], exons2[m][0]) - min(exons1[n][1], exons2[m][1])) + 1
                    for n in range(len(exons1))
                    for m in range(len(exons2))
                    if exons1[n][0] <= exons2[m][1] and exons1[n][1] >= exons2[m][0]
                ])
                
                # Calculate the percentage of each sequence being overlapped by the other
                len1 = sum([ end - start + 1 for start, end in exons1 ])
                len2 = sum([ end - start + 1 for start, end in exons2 ])
                
                pct1 = exonOverlap / len1
                pct2 = exonOverlap / len2
                
                # See if this meets any coverage pct cutoffs
                if len1 <= len2:
                    shorterPct = pct1
                    longerPct = pct2
                else:
                    shorterPct = pct2
                    longerPct = pct1
                
                if shorterPct >= shorterCovPct and longerPct >= longerCovPct:
                    # If this met the cutoff, add an edge between these nodes in the graph
                    idsGraph.add_edge(id1, id2)
        
        # Identify connected IDs which form clusters
        binList = []
        for connectedIDs in nx.connected_components(idsGraph):
            start = 1 # start is not needed anymore
            end = 10 # end is not needed anymore
            exonsDict = { seqID : [] for seqID in connectedIDs } # exons are not needed anymore
            
            # Create a bin for this split cluster
            newBin = Bin(bin.contig, start, end)
            newBin.union(connectedIDs, exonsDict)
            
            binList.append(newBin)
        
        return binList
