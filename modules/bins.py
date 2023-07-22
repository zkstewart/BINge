from intervaltree import IntervalTree, Interval
import networkx as nx

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
    
    def __repr__(self):
        return (f"<Bin object;contig='{self.contig}';start={self.start};" +
                f"end={self.end};num_ids={len(self.ids)}"
        )
    
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
    
    def __len__(self):
        return len(self.bins)
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinCollection object;num_bins={0}>".format(
            len(self)
        )

class BinSplitter:
    '''
    The BinSplitter Class provides a simple tool for clustering sequences
    within a bin (which essentially is splitting the bin up)
    based on whether they have overlaps over any probable exon regions.
    It isn't intended to be effective for most applications, but
    it may prove useful for BINge to address a particular need that CD-HIT
    does not (i.e., simply checking for whether two sequences have any shared
    regions or not!).
    
    Attributes:
        bin (REQUIRED) -- a Bin object as defined in this module.
        shorter_cov_pct (OPTIONAL) -- a float from 0.0 -> 1.0 indicating how much
                                      the shorter sequence must cover a longer sequence
                                      for them to be considered as sharing sequence
        longer_cov_pct (OPTIONAL) -- a float from 0.0 -> 1.0 indicating how much
                                     the longer sequence must be covered by a shorter
                                     sequence for them to be considered as sharing sequence
    '''
    def __init__(self, bin):
        # Set attributes based on object init arguments
        self.bin = bin
        
        # Set default attributes
        self.shorter_cov_pct = 0.25
        self.longer_cov_pct = 0.01
        
        # Results storage values
        self.resultClusters = None
    
    @property
    def shorter_cov_pct(self):
        return self._shorter_cov_pct
    
    @shorter_cov_pct.setter
    def shorter_cov_pct(self, num):
        '''
        Controls sequence clustering based on overlap percentage for the
        shorter sequence.
        
        Parameters:
            num -- a float from 0.0 -> 1.0 indicating how much
                   the shorter sequence must cover a longer sequence
                   for them to be considered as sharing sequence
        '''
        assert isinstance(num, float) or isinstance(num, int)
        assert 0.0 <= num <= 1.0, \
            "shorter_cov_pct must be in the range bounded by 0 and 1"
        self._shorter_cov_pct = num
    
    @property
    def longer_cov_pct(self):
        return self._longer_cov_pct
    
    @longer_cov_pct.setter
    def longer_cov_pct(self, num):
        '''
        Controls sequence clustering based on overlap percentage for the
        longer sequence.
        
        Parameters:
            num -- a float from 0.0 -> 1.0 indicating how much
                   the longer sequence must be covered by a shorter sequence
                   for them to be considered as sharing sequence
        '''
        assert isinstance(num, float) or isinstance(num, int)
        assert 0.0 <= num <= 1.0, \
            "longer_cov_pct must be in the range bounded by 0 and 1"
        self._longer_cov_pct = num
    
    def cluster(self):
        '''
        The main method of the ExonOverlap Class. This function will automatically cluster
        sequences based on whether they have some degree of exon overlap.
        
        Sets:
            self.resultClusters -- will be set to a dictionary with structure like:
                                   {
                                       0: [
                                           seqID1, seqID2, ...
                                       ],
                                       1: [
                                           seqID#, seqID#, ...
                                       ],
                                       ...
                                   }
        '''
        ids = list(self.bin.ids) # ensure we always use this in consistent ordering
        
        # Model links between sequences as a graph structure
        idsGraph = nx.Graph()
        idsGraph.add_nodes_from(ids)
        
        # Figure out which sequences have links through shared exons
        for i in range(len(ids)-1):
            for x in range(i+1, len(ids)):
                id1, id2 = ids[i], ids[x]
                exons1, exons2 = self.bin.exons[id1], self.bin.exons[id2]
                
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
                
                if shorterPct >= self.shorter_cov_pct and longerPct >= self.longer_cov_pct:
                    # If this met the cutoff, add an edge between these nodes in the graph
                    idsGraph.add_edge(id1, id2)
        
        # Identify connected IDs which form clusters
        ongoingCount = 0
        clusterDict = {}
        for connectedIDs in nx.connected_components(idsGraph):
            connectedIDs = list(connectedIDs)
            clusterDict[ongoingCount] = connectedIDs
            ongoingCount += 1
        
        # Set value
        self.resultClusters = clusterDict
