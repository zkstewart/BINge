from intervaltree import IntervalTree, Interval
import networkx as nx
from hashlib import sha256

class Bin:
    '''
    Make sure to use this bin with 1-based indexing for start and end values,
    in the same way as a GFF3 would.
    '''
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start
        self.end = end
        
        self.ids = set() # set of sequence IDs
    
    @property
    def contig(self):
        return self._contig
    
    @contig.setter
    def contig(self, value):
        assert isinstance(value, str)
        assert value != ""
        
        self._contig = value
    
    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self, value):
        assert isinstance(value, int)
        assert value > 0, \
            "A Bin must start at a value >= 1 (it behaves 1-based for indexing)"
        self._start = value
    
    @property
    def end(self):
        return self._end
    
    @end.setter
    def end(self, value):
        assert isinstance(value, int)
        assert value > 0, \
            "A Bin must end at a value >= 1 (it behaves 1-based for indexing)"
        self._end = value
    
    def add(self, idValue):
        assert isinstance(idValue, str)
        self.ids.add(idValue)
    
    def union(self, ids):
        assert isinstance(ids, list) or isinstance(ids, set)
        self.ids = self.ids.union(ids)
    
    def merge(self, otherBin):
        assert self.contig == otherBin.contig, \
            "Cannot merge bins on different contigs!"
        
        self.start = min(self.start, otherBin.start)
        self.end = max(self.end, otherBin.end)
        self.union(otherBin.ids)
    
    def is_overlapping(self, otherBin, shorterCovPct=0.25, longerCovPct=0.10):
        '''
        Compares this bin to another bin to see if they overlap each other according
        to shorter sequence and longer sequence coverage cut-offs. Think of this
        similarly to how CD-HIT has the -aS and -aL values.
        
        Parameters:
            otherBin -- a different Bin object to compare.
            shorterCovPct / longerCovPct -- a float value in the range of 0 -> 1 which
                                            controls whether any sequences in the two bins
                                            are deemed to overlap or not. Overlap occurs
                                            when both values are satisfied (AND), not when
                                            a single value is satisfied (OR).
        '''
        for num in [shorterCovPct, longerCovPct]:
            assert isinstance(num, float) or isinstance(num, int)
            assert 0.0 <= num <= 1.0, \
                "shorterCovPct and longerCovPct must be in the range bounded by 0 and 1"
        
        if not (self.start <= otherBin.end and self.end >= otherBin.start):
            return False
        else:
            # Calculate the percentage of overlap between the two bins
            overlap = min(self.end, otherBin.end) - max(self.start, otherBin.start) + 1
            assert overlap > 0, "Logic error in overlap calculation"
            
            len1 = self.end - self.start + 1
            len2 = otherBin.end - otherBin.start + 1
            
            pct1 = overlap / len1
            pct2 = overlap / len2
            
            # See if this meets any coverage pct cutoffs
            if len1 <= len2:
                shorterPct = pct1
                longerPct = pct2
            else:
                shorterPct = pct2
                longerPct = pct1
            
            if shorterPct >= shorterCovPct and longerPct >= longerCovPct:
                return True
            else:
                return False
        
    def __repr__(self):
        return (f"<Bin object;contig='{self.contig}';start={self.start};" +
                f"end={self.end};num_ids={len(self.ids)}"
        )
    
    def __hash__(self):
        return hash((self.contig, self.start, self.end,
                     tuple(self.ids)))
    
    def sha256(self):
        h = sha256()
        for s in (
                self.contig.encode(), str(self.start).encode(), str(self.end).encode(),
                str(self.ids).encode()
        ):
            h.update(s)

        return h.hexdigest()
    
    def __eq__(self, other):
        if isinstance(other, Bin):
            return self.sha256() == other.sha256()
        return False

class BinCollection:
    '''
    Encapsulates an IntervalTree allowing easy addition and finding of Bin objects.
    Indexes by exon not CDS.
    '''
    def __init__(self):
        self.bins = {}
    
    def add(self, bin):
        self.bins.setdefault(bin.contig, IntervalTree())
        self.bins[bin.contig][bin.start:bin.end+1] = bin
    
    def find(self, contig, start, end):
        '''
        Start and end should be provided as 1-based and inclusive values.
        For example, searching for 100,100 will find overlaps at that exact
        position.
        '''
        if not contig in self.bins:
            return []
        else:
            bins = [ b.data for b in self.bins[contig][start:end+1] ]
        return bins
    
    def delete(self, bin):
        if not bin.contig in self.bins:
            raise ValueError(f"This BinCollection does not contain any bins on contig '{bin.contig}'")
        else:
            self.bins[bin.contig].remove(Interval(bin.start, bin.end+1, bin))
    
    def replace(self, oldBin, newBin):
        self.delete(oldBin)
        self.add(newBin)
    
    def merge(self, otherBinCollection):
        '''
        Merges the otherBinCollection into this one, adding all its Bins into this. This
        results in changes to this object.
        
        This will NOT merge overlapping bins in any way. It just purely adds the Bins of
        the other collection into this one.
        '''
        for contig, tree in otherBinCollection.bins.items():
            for bin in tree:
                self.add(bin.data)
    
    def __len__(self):
        return sum([ len(tree) for contig, tree in self.bins.items() ])
    
    def __iter__(self):
        for contig, tree in self.bins.items():
            for bin in tree:
                yield bin
    
    def __repr__(self):
        return "<BinCollection object;num_bins={0}>".format(
            len(self)
        )

class BinBundle:
    '''
    Encapsulates a simple list of Bin objects for storage and manipulation when the searching
    of an IntervalTree (as provided by BinCollection) is not needed.
    
    The naming is arbitrary. A bundle is less organised than a collection I suppose?
    '''
    def __init__(self):
        self.bins = []
    
    def add(self, bin):
        self.bins.append(bin)
    
    def merge(self, otherBinCollection, otherType="BinCollection"):
        '''
        Merges the otherBinCollection into this one, adding all its Bins into this. This
        results in changes to this object. Can handle BinBundle or BinCollection
        objects.
        '''
        assert otherType in ["BinCollection", "BinBundle"], \
            "otherType must be either 'BinCollection' or 'BinBundle'"
        
        if otherType == "BinCollection":
            for bin in otherBinCollection:
                self.add(bin.data)
        else:
            for bin in otherBinCollection:
                self.add(bin)
    
    def link_bins(self, VOTE_THRESHOLD = 0.5):
        '''
        This will attempt to merge bins that are "equivalent" across the genomes using
        a graph-based approach to modelling the relationships between the bins. The
        result is a BinBundle that unifies all input data and is appropriate for
        sequence clustering.
        
        Note: this function is not guaranteed to produce a stable output i.e., if
        you feed its result back in, you may get different results. It should
        eventually converge on a stable solution, but there might theoretically
        be an edge case where that never happens.
        
        Parameters:
            VOTE_THRESHOLD -- a float greater than 0, and less than or equal to 1.0.
                              This represents a ratio 
        Returns:
            linkedBinBundle -- a NEW BinBundle object where bins from the
                               current object have been merged where deemed
                               appropriate.
        '''
        assert 0 < VOTE_THRESHOLD <= 1.0, \
            "VOTE_THRESHOLD must be a value greater than 0, and less than or equal to 1"
        
        # Format bins into a dictionary
        "It's not ideal for memory use but I need instant lookup and consistent ordering"
        binDict = {}
        binCounter = 0
        for bin in self:
            binDict[binCounter] = bin
            binCounter += 1
        
        # Figure out which bins have links through shared sequences
        idLinks = {}
        for binIndex, bin in binDict.items():
            for seqID in bin.ids:
                idLinks.setdefault(seqID, set())
                idLinks[seqID].add(binIndex)
        
        # Model links between bins as a graph structure
        binGraph = nx.Graph()
        binGraph.add_nodes_from(range(0, len(binDict)))
        
        for binIndex, bin in binDict.items():
            # Perform tallying and identify edges in the bin merging graph
            votes = {}
            for seqID in bin.ids:
                for index in idLinks[seqID]:
                    if index != binIndex:
                        votes.setdefault(index, 0)
                        votes[index] += 1
            toLink = [
                index
                for index, count in votes.items()
                if count / len(bin.ids) >= VOTE_THRESHOLD
            ]
            
            # Add edges between any relevant bin nodes
            if len(toLink) > 0:
                for index in toLink:
                    binGraph.add_edge(binIndex, index)
        
        # Merge bins on the basis of link identification
        "At this point, the bins will lose any meaning they have in their .contig, .start, etc"
        linkedBinBundle = BinBundle()
        for connectedBins in nx.connected_components(binGraph):
            connectedBins = list(connectedBins)
            
            newBin = binDict[connectedBins[0]]
            for index in connectedBins[1:]:
                newBin.union(binDict[index].ids)
            linkedBinBundle.add(newBin)
        
        return linkedBinBundle
    
    @staticmethod
    def create_from_collection(binCollection):
        '''
        Initializes a BinBundle object from a BinCollection object.
        
        Parameters:
            binCollection -- a BinCollection object.
        Returns:
            binBundle -- a BinBundle object.
        '''
        newBundle = BinBundle()
        for bin in binCollection:
            newBundle.add(bin.data)
        return newBundle
    
    @staticmethod
    def create_from_multiple_collections(binCollectionList):
        '''
        Initializes a BinBundle object from one or more BinCollection objects.
        The input value is expected to be a list of BinCollection objects.
        The result is a single BinBundle object with all the bins squished
        down into it.
        
        Parameters:
            binCollectionList -- a list of BinCollection objects.
        Returns:
            binBundle -- a BinBundle object.
        '''
        assert isinstance(binCollectionList, list), \
            "binCollectionList must be a list of BinCollection objects"
        
        newBundle = BinBundle()
        for binCollection in binCollectionList:
            for bin in binCollection:
                newBundle.add(bin.data)
        return newBundle
        
        
    def __len__(self):
        return len(self.bins)
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinBundle object;num_bins={0}>".format(
            len(self)
        )
