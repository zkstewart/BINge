import os, sys
from multiprocessing import Process, Pipe, Queue

from .gff3_handling import GFF3, iterate_gmap_gff3
from .fasta_handling import load_sequence_length_index
from .bins import Bin, BinCollection

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
    
    # If it overlaps only one exon, merge it in normally
    elif len(binOverlap) == 1:
        overlappingBin = binOverlap[0]

        # ... and add the new bin in its place
        newBin.merge(overlappingBin)
        binCollection.delete(overlappingBin)
        binCollection.add(newBin)
    
    # It it has multiple overlaps, add its ID in without a true merge
    else:
        for overlappingBin in binOverlap:
            overlappingBin.union(newBin.ids)

def find_overlapping_bins(binCollection, binQuery):
    '''
    Assistant function to help with taking the bins found to overlap a query region
    via the IntervalTree underlying a BinCollection, and then filtering those bins
    to only those that share a substantial proportion of overlap.
    
    Parameters:
        binCollection -- a BinCollection() object to use the .find() method on.
        binQuery -- a Bin() which we want to find overlapping bins from the binCollection
                    from.
    Returns:
        binOverlap -- a list of Bin() objects that overlap the binQuery, and share a
                      substantial proportion of overlap.
    '''
    foundBins = binCollection.find(binQuery.contig, binQuery.start, binQuery.end)
    
    binOverlap = []
    for bin in foundBins:
        if binQuery.is_overlapping(bin):
            binOverlap.append(bin)
    return binOverlap

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
        isMicrobial -- a boolean indicating whether the genomes are microbial or not which,
                       in turn, determines whether we will parse mRNA features (False) or
                       gene features (True).
    '''
    def task(self, gff3File, isMicrobial=False):
        binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if gff3File != None:
            gff3Obj = GFF3(gff3File, strict_parse=False)
            for geneFeature in gff3Obj.types["gene"]:
                
                # Create a bin for each exon feature
                exonBins = []
                if not isMicrobial:
                    try:
                        for mrnaFeature in geneFeature.mRNA:
                            for exonFeature in mrnaFeature.exon:
                                exonBin = Bin(exonFeature.contig, exonFeature.start, exonFeature.end)
                                exonBin.add(mrnaFeature.ID)
                                exonBins.append(exonBin)
                    except:
                        "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                        continue
                else:
                    try:
                        exonBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
                        exonBin.add(geneFeature.ID)
                        exonBins.append(exonBin)
                    except:
                        "This exception occurs if a gene feature has unusual children types e.g., rRNA"
                        continue
                
                # Iteratively handle exon bins
                for exonBin in exonBins:
                    binOverlap = find_overlapping_bins(binCollection, exonBin)
                    add_bin_to_collection(binCollection, binOverlap, exonBin)
        
        return binCollection

class GmapBinProcess(ReturningProcess):
    '''
    Bins GMAP alignments into each genome's BinCollection object.
    
    Parameters:
        gmapFile -- a string indicating the location of a GMAP GFF3 for parsing.
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        indexFileName -- a string indicating the location of a pickled dictionary
                         linking sequence IDs (key) to lengths (value).
        minIdentity -- (optional) a float fraction indicating what identity value is
                       minimally required for us to use a GMAP alignment; default=0.95.
    '''
    def task(self, gmapFiles, binCollection, indexFileName, minIdentity=0.95):
        # Behavioural parameters (static for now, may change later)
        "These statics are for filtering GMAP alignments that are poor quality"
        OKAY_COVERAGE = 96.5
        OKAY_INDEL_PROPORTION = 0.01
        ##
        "These statics prevent a gene mapping to multiple locations when it has a clear best"
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        ##
        "This static allows coverage to be lower if the alignment is close to the contig's edge"
        LENIENT_BOUNDARY_LENGTH = 1000
        
        # Load in the lengths index
        seqLenDict = load_sequence_length_index(indexFileName)
        
        # Iterate through GMAP files
        for gmapFile in gmapFiles:
            # Hold onto the GMAP alignments for each gene
            """GMAP is guaranteed to return alignments for each gene together, but not ordered
            by quality. We'll sort them by quality and then process them iteratively"""
            thisBlockID, thisBlockData = None, None
            for dataDict in iterate_gmap_gff3(gmapFile):
                # Handle first iteration
                if thisBlockID == None:
                    thisBlockID = dataDict["Name"]
                    thisBlockData = [dataDict]
                    continue
                
                # Store data if this is from the same sequence
                if dataDict["Name"] == thisBlockID:
                    thisBlockData.append(dataDict)
                    continue
                
                # Otherwise, sort the data dicts by their quality
                "Higher coverage -> higher identity -> lower indels"
                thisBlockData.sort(key = lambda x: (-x["coverage"], -x["identity"], x["indels"]))
                
                # Iterate through data dicts
                bestCov, bestIdent = None, None
                for blockDataDict in thisBlockData:
                    # Get alignment statistics
                    coverage, identity = blockDataDict["coverage"], blockDataDict["identity"]
                    exonLength = sum([
                        end - start + 1
                        for start, end in blockDataDict["exons"]
                    ])
                    indelProportion = blockDataDict["indels"] / exonLength
                    
                    # Apply leniency if the alignment is close to the contig's boundaries
                    """This helps to deal with fragmented genes aligning well to a contig's
                    edge, but the full length model not being binned because its coverage is
                    too low."""
                    beLenient = False
                    if coverage < OKAY_COVERAGE:
                        contigLength = seqLenDict[blockDataDict["contig"]]
                        alignmentStart = min(min(blockDataDict["exons"]))
                        alignmentEnd = max(max(blockDataDict["exons"]))
                        
                        if (alignmentStart < LENIENT_BOUNDARY_LENGTH) \
                            or (alignmentEnd + LENIENT_BOUNDARY_LENGTH) > contigLength:
                            beLenient = True
                    
                    # Skip processing if the alignment sucks
                    isGoodAlignment = (True if beLenient else coverage >= OKAY_COVERAGE) \
                                    and identity >= minIdentity \
                                    and indelProportion <= OKAY_INDEL_PROPORTION
                    if not isGoodAlignment:
                        continue
                    
                    # See if this path should be skipped because another better one was processed
                    if bestCov != None:
                        if (coverage + ALLOWED_COV_DIFF) < bestCov \
                            or (identity + ALLOWED_IDENT_DIFF) < bestIdent:
                            continue
                    else:
                        bestCov, bestIdent = coverage, identity # set if this is the 'best'
                    
                    # Iteratively handle exon features
                    for exonStart, exonEnd in blockDataDict["exons"]:
                        # Create a bin for each exon feature
                        exonBin = Bin(blockDataDict["contig"], exonStart, exonEnd)
                        exonBin.add(blockDataDict["Name"])
                        
                        # See if this overlaps an existing bin
                        binOverlap = find_overlapping_bins(binCollection, exonBin)
                        
                        # Add a new bin, or merge any bins as appropriate
                        add_bin_to_collection(binCollection, binOverlap, exonBin)
                
                # Reset for next iteration
                thisBlockID, thisBlockData = dataDict["Name"], [dataDict]
        
        return binCollection

class GraphPruneProcess(ReturningProcess):
    '''
    Allows for BinGraph.prune() to be run in parallel for multiple genomes.
    
    Parameters:
        binGraph -- a BinGraph object to call .prune() on.
        WEIGHT_CUTOFF -- (optional) a float indicating the minimum weight proportion
                         required for an edge to be kept; default=0.5.
        CUT_CUTOFF -- (optional) a float indicating the minimum weight proportion
                       required for a sequence ID to not be eliminated if it shows
                       up in cut edges; default=0.5.
    Returns:
        binGraph -- the same BinGraph object, but with its .prune() method called.
    '''
    def task(self, binGraph, WEIGHT_CUTOFF=0.5):
        chimers = binGraph.find_chimers_by_pruning(WEIGHT_CUTOFF)
        return chimers
