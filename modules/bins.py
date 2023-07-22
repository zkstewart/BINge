from intervaltree import IntervalTree, Interval

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
                exonFeature.coords
                for mrnaFeature in gff3Feature.mRNA
                for exonFeature in mrnaFeature.exon
            ])
            tree.merge_overlaps()
            exons = sorted([
                [interval.begin, interval.end]
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
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinCollection object;num_bins={0}>".format(
            len(self.bins)
        )
