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
        self.exons = {} # dictionary pairing set IDs to exon lists
    
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
    
    def add(self, idValue, exonList):
        assert isinstance(idValue, str)
        self.ids.add(idValue)
        self.exons[idValue] = exonList
    
    def union(self, ids, exons):
        assert isinstance(ids, list) or isinstance(ids, set)
        assert isinstance(exons, dict)
        self.ids = self.ids.union(ids)
        self.exons.update(exons)
    
    def merge(self, otherBin):
        assert self.contig == otherBin.contig, \
            "Cannot merge bins on different contigs!"
        
        self.start = min(self.start, otherBin.start)
        self.end = max(self.end, otherBin.end)
        self.union(otherBin.ids, otherBin.exons)
    
    def is_overlapping(self, otherBin, shorterCovPct=0.25, longerCovPct=0.01):
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
        
        for thisExon in self.exons.values():
            for otherExon in otherBin.exons.values():
                exonOverlap = sum([
                    abs(max(thisExon[n][0], otherExon[m][0]) - min(thisExon[n][1], otherExon[m][1])) + 1
                    for n in range(len(thisExon))
                    for m in range(len(otherExon))
                    if thisExon[n][0] <= otherExon[m][1] and thisExon[n][1] >= otherExon[m][0]
                ])
                
                # Calculate the percentage of each sequence being overlapped by the other
                len1 = sum([ end - start + 1 for start, end in thisExon ])
                len2 = sum([ end - start + 1 for start, end in otherExon ])
                
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
                    return True
        return False
    
    def __repr__(self):
        return (f"<Bin object;contig='{self.contig}';start={self.start};" +
                f"end={self.end};num_ids={len(self.ids)}"
        )
    
    def __hash__(self):
        return hash((self.contig, self.start, self.end,
                     tuple(self.ids), hash(tuple(map(tuple, self.exons)))))
    
    def sha256(self):
        h = sha256()
        for s in (
                self.contig.encode(), str(self.start).encode(), str(self.end).encode(),
                str(self.ids).encode(), str(self.exons).encode()
        ):
            h.update(s)

        return h.hexdigest()
    
    def __eq__(self, other):
        if isinstance(other, Bin):
            return self.sha256() == other.sha256()
        return False
    
    @staticmethod
    def format_exons_from_gff3_feature(gff3Feature):
        '''
        Helper function for getting exons from a GFF3 feature when creating
        a Bin object. For mRNAs, this is simple - we just get the .exon coordinates.
        For genes, we use the IntervalTree class to merge overlapping coordinate
        ranges.
        
        Parameters:
            gff3Feature -- a GFF3.Feature for a gene or mRNA object
        '''
        assert gff3Feature.type in ["gene", "mRNA"], \
            "Can only format exons from a GFF3 feature that is a gene or mRNA!"
        
        if gff3Feature.type == "gene":
            assert hasattr(gff3Feature, "mRNA"), \
                "Gene GFF3 feature does not have mRNA attribute; cannot format exons..."
            
            tree = IntervalTree.from_tuples([
                [ exonFeature.start, exonFeature.end+1 ] # prevent errors for 1-bp exons
                for mrnaFeature in gff3Feature.mRNA
                for exonFeature in mrnaFeature.exon
            ])
            tree.merge_overlaps()
            exons = sorted([
                [interval.begin, interval.end-1]
                for interval in tree
            ])
        else:
            exons = sorted([
                exonFeature.coords
                for exonFeature in gff3Feature.exon
            ])

        return exons

class BinCollection:
    '''
    Encapsulates an IntervalTree allowing easy addition and finding of Bin objects.
    Indexes by exon not CDS.
    '''
    def __init__(self):
        self.bins = IntervalTree()
    
    def add(self, bin):
        self.bins[bin.start:bin.end+1] = bin
    
    def find(self, contig, start, end):
        '''
        Start and end should be provided as 1-based and inclusive values.
        For example, searching for 100,100 will find overlaps at that exact
        position.
        '''
        bins = [ b.data for b in self.bins[start:end+1] if b.data.contig == contig ]
        return bins
    
    def delete(self, bin):
        self.bins.remove(Interval(bin.start, bin.end+1, bin))
    
    def replace(self, oldBin, newBin):
        self.bins.remove(Interval(oldBin.start, oldBin.end+1, oldBin))
        self.add(newBin)
    
    def merge(self, otherBinCollection):
        '''
        Merges the otherBinCollection into this one, adding all its Bins into this. This
        results in changes to this object.
        
        This will NOT merge overlapping bins in any way. It just purely adds the Bins of
        the other collection into this one.
        '''
        for bin in otherBinCollection:
            self.add(bin.data)
    
    def flatten(self, otherBinCollection):
        '''
        Flattens the otherBinCollection into this one, adding all its Bins into this. This
        results in changes to this object.
        
        As opposed to .merge(), this WILL merge overlapping bins. It's a more costly operation
        so it should only be done if necessary.
        '''
        for bin in otherBinCollection:
            ovlBins = self.find(bin.data.contig, bin.data.start, bin.data.end)
            if len(ovlBins) == 0:
                self.add(bin.data)
            else:
                newBin = ovlBins[0]
                for ovlBin in ovlBins[1:]:
                    newBin.merge(ovlBin)
                newBin.merge(bin.data)
    
    def __len__(self):
        return len(self.bins)
    
    def __iter__(self):
        return iter(self.bins)
    
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
                newBin.union(binDict[index].ids, binDict[index].exons)
            linkedBinBundle.add(newBin)
        
        return linkedBinBundle
    
    def __len__(self):
        return len(self.bins)
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinBundle object;num_bins={0}>".format(
            len(self)
        )
