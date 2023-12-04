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
    
    def link_bins(self, VOTE_THRESHOLD = 0.5):
        '''
        This will attempt to merge bins that are "equivalent" across the genomes using
        a graph-based approach to modelling the relationships between the bins. The
        result is a BinCollection that unifies all input data and is appropriate for
        sequence clustering.
        
        After performing this process, you should expect that the bin attributes
        relating to sequence coordinates i.e., contig, start, and end, will no
        longer be accurate. *Don't try to use them afterwards!*
        
        Note: this function is not guaranteed to produce a stable output i.e., if
        you feed its result back in, you may get different results. It should
        eventually converge on a stable solution, but there might theoretically
        be an edge case where that never happens.
        
        Parameters:
            VOTE_THRESHOLD -- a float greater than 0, and less than or equal to 1.0.
                              This represents a ratio 
        Returns:
            linkedBinCollection -- a NEW BinCollection object where bins from the
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
            binDict[binCounter] = bin.data
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
        linkedBinCollection = BinCollection()
        for connectedBins in nx.connected_components(binGraph):
            connectedBins = list(connectedBins)
            
            newBin = binDict[connectedBins[0]]
            for index in connectedBins[1:]:
                newBin.union(binDict[index].ids, binDict[index].exons)
            linkedBinCollection.add(newBin)
        
        return linkedBinCollection
    
    def fix_fragments(self, multiOverlap, VOTE_THRESHOLD = 0.5):
        '''
        Attempt to merge bins that are likely to be fragmentary.
        
        It does this using 'multiOverlap', a list of GMAP alignment Features which
        overlap more than one Bin. It is assumed that these Bins have been previously
        seeded from the genome's annotation, and hence this function will help us to
        identify situations where the genome annotation has flaws.
        
        The method of performing this is similar to that seen in bin_self_linker()
        where it's leveraging the graph-based approach to model bins which should
        merge. Merging decisions are made based on comparing the 'weight' of an edge
        (i.e., how many GMAP alignments support the bin's merging) to the number of
        IDs which do not support the merging (i.e., how many GMAP alignments only
        overlapped a single bin). The vote threshold decides whether we should merge
        or not.
        
        Note: it's unsure if this function will produce a stable output i.e., it's
        possible that, if you feed its result back in, you may get different results.
        Theoretically, the first pass of this function should capture just about
        everything so it probably doesn't need to be run iteratively.
        
        Parameters:
            multiOverlap -- a list containing GFF3 Features which were found to
                            overlap more than one Bin in the binCollection.
            VOTE_THRESHOLD -- a float value in the range of 0 < VOTE_THRESHOLD <= 1.0
                              determining what percentage (e.g., 0.5 = 50%) of GMAP
                              alignments should support bin merging via overlap before
                              we actually merge the bins together.
        Returns:
            linkedBinCollection -- a NEW BinCollection object where bins from the
                                   current object have been merged where deemed
                                   appropriate.
        '''
        assert 0 < VOTE_THRESHOLD <= 1.0, \
            "VOTE_THRESHOLD must be a value greater than 0, and less than or equal to 1"
        
        # Format bins into a dictionary
        binDict = {}
        numBins = 0
        for bin in self:
            binDict[bin.data.sha256()] = bin.data
            numBins += 1
        assert len(binDict) == numBins, \
            "Hash collision occurred! I can't handle this, and you should buy a lottery ticket..."
        
        # Model links between bins as a graph structure
        binGraph = nx.Graph()
        binGraph.add_nodes_from(binDict.keys())
        
        # Process bins that are multi-overlapped
        for overlappingFeature in multiOverlap:
            # Find the bins which are overlapped
            binOverlap = self.find(overlappingFeature.contig, overlappingFeature.start, overlappingFeature.end)
            assert len(binOverlap) > 1, "Found a non-multi-overlapper somehow?"
            binOverlap = [ bin.sha256() for bin in binOverlap ]
            
            # Add edges between overlapped bins
            for i in range(0, len(binOverlap)-1):
                for x in range(i+1, len(binOverlap)):
                    baseID = overlappingFeature.ID.rsplit(".", maxsplit=1)[0]
                    exons = Bin.format_exons_from_gff3_feature(overlappingFeature)
                    
                    # Add new edge
                    if not binGraph.has_edge(binOverlap[i], binOverlap[x]):
                        binGraph.add_edge(binOverlap[i], binOverlap[x],
                                        weight=1, ids=[baseID],
                                        exons={baseID: exons})
                    
                    # Increase weight of existing edge
                    else:
                        binGraph[binOverlap[i]][binOverlap[x]]["weight"] += 1
                        binGraph[binOverlap[i]][binOverlap[x]]["ids"].append(baseID)
                        binGraph[binOverlap[i]][binOverlap[x]]["exons"][baseID] = exons
        
        # Create a new network with edges set if they exceed the unlinked number of IDs
        finalGraph = nx.Graph()
        finalGraph.add_nodes_from(binDict.keys())
        
        for connectedBins in nx.connected_components(binGraph):
            connectedBins = list(connectedBins)
            
            # Check that this bin is fully connected
            "Non-fully connected graphs at this stage are dodgy"
            allConnected = True
            for i in range(0, len(connectedBins)-1):
                for x in range(i+1, len(connectedBins)):
                    binHash1 = connectedBins[i]
                    binHash2 = connectedBins[x]
                    try:
                        binGraph[binHash1][binHash2]
                    except:
                        allConnected = False
                        break # small speed up, won't fully exit loop
            if allConnected == False:
                continue
            
            # Check weights for voting and making connections between nodes
            for i in range(0, len(connectedBins)-1):
                for x in range(i+1, len(connectedBins)):
                    binHash1 = connectedBins[i]
                    binHash2 = connectedBins[x]
                    edgeWeight = binGraph[binHash1][binHash2]["weight"]
                    distanceWeight = sum([
                        len(binDict[binHash1].ids), len(binDict[binHash2].ids)
                    ])
                    
                    mergeVote = edgeWeight / (edgeWeight + distanceWeight)
                    if mergeVote >= VOTE_THRESHOLD:
                        finalGraph.add_edge(binHash1, binHash2,
                                            ids = binGraph[binHash1][binHash2]["ids"],
                                            exons = binGraph[binHash1][binHash2]["exons"])
        
        # Merge bins on the basis of link identification
        linkedBinCollection = BinCollection()
        for connectedBins in nx.connected_components(finalGraph):
            connectedBins = list(connectedBins)
            
            newBin = binDict[connectedBins[0]]
            binHash1 = connectedBins[0]
            for binHash2 in connectedBins[1:]:
                newBin2 = Bin(newBin.contig, binDict[binHash2].start, binDict[binHash2].end)
                newBin2.ids = set(finalGraph[binHash1][binHash2]["ids"]).union(binDict[binHash2].ids)
                
                newBin2.exons = finalGraph[binHash1][binHash2]["exons"]
                newBin2.exons.update(binDict[binHash2].exons)
                
                newBin.merge(newBin2)
            linkedBinCollection.add(newBin)
        
        return linkedBinCollection
    
    def fix_chimeras(self, VOTE_THRESHOLD = 0.5):
        '''
        Attempt to fragment bins that are likely to be chimeric. This is different
        than what BinSplitter does, as it's not merely separating bin members into
        2 or more individual bins. Instead, it's cutting one or more lines through
        a bin, and ...
        
        I think this function is too hard to implement here. It's really deserving
        of being its own fully fledged program.
        
        Parameters:
            VOTE_THRESHOLD -- a float value in the range of 0 < VOTE_THRESHOLD <= 1.0
                              determining what percentage (e.g., 0.5 = 50%) of GMAP
                              alignments should support bin merging via overlap before
                              we actually merge the bins together.
        Returns:
            linkedBinCollection -- a NEW BinCollection object where bins from the
                                   current object have been merged where deemed
                                   appropriate.
        '''
        assert 0 < VOTE_THRESHOLD <= 1.0, \
            "VOTE_THRESHOLD must be a value greater than 0, and less than or equal to 1"
        
        raise NotImplementedError()
    
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
