import os, sys
from multiprocessing import Process, Pipe, Queue

from .gff3_handling import GFF3, iterate_through_gff3
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
                                and identity >= minIdentity \
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
                coordinateOverlaps = binCollection.find(feature.contig, feature.start, feature.end)
                binOverlap = []
                for bin in coordinateOverlaps:
                    shouldMerge = newBin.is_overlapping(bin)
                    if shouldMerge:
                        binOverlap.append(bin)
                
                # If it does not overlap an existing bin ...
                if len(binOverlap) == 0:
                    # ... add the novel bin into the collection
                    binCollection.add(newBin)
                
                # Exclude any transcripts that overlap multiple bins
                elif len(binOverlap) > 1:
                    """We expect these occurrences to be chimeric transcripts which should
                    be filtered. Otherwise, fragmented transcripts being loaded in may
                    result in a single loci having multiple bins over it if the full length
                    one doesn't get processed first. We'll need to fix those later by
                    assessing these multi-overlaps."""
                    multiOverlaps.append(feature)
                
                # If it overlaps exactly 1 bin ...
                else:
                    # ... merge the bins together
                    add_bin_to_collection(binCollection, binOverlap, newBin)
        
        return binCollection, multiOverlaps

class FragmentFixProcess(ReturningProcess):
    '''
    Allows for fragment fixing to be run in parallel across multiple separate
    BinCollection() objects.
    
    Parameters:
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        multiOverlaps -- a list containing GFF3 Feature objects indicating
                         GMAP alignments that overlapped more than one bin.
    '''
    def task(self, binCollection, multiOverlaps):
        binCollection.fix_fragments(multiOverlaps)
        return binCollection
