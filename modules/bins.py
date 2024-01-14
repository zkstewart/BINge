from intervaltree import IntervalTree, Interval
import networkx as nx
from hashlib import sha256

class Bin:
    '''
    Make sure to use this bin with 1-based indexing for start and end values,
    in the same way as a GFF3 would.
    
    Parameters:
        contig -- a string indicating the contig name that this bin is on.
        start -- an integer indicating the start position of the bin.
        end -- an integer indicating the end position of the bin.
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
    
    def is_overlapping(self, otherBin, shorterCovPct=0.20, longerCovPct=0.10):
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
            
            if shorterPct >= shorterCovPct or longerPct >= longerCovPct:
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

class BinGraph:
    '''
    Represents a graph where nodes are Bin objects and edges are the sequence alignments
    that join them.
    '''
    def __init__(self, binCollection):
        self.graph = nx.Graph()
        self.graph.add_nodes_from([ interval.data for interval in binCollection ])
        self._link_bins_via_collection(binCollection)
        self.eliminations = set()
    
    def _link_bins_via_collection(self, binCollection):
        # Locate all bins that each sequence ID occurs in
        idLinks = {}
        for interval in binCollection:
            bin = interval.data
            for seqID in bin.ids:
                idLinks.setdefault(seqID, [])
                idLinks[seqID].append(bin)
        
        # Go through and link bins together
        for seqID, binList in idLinks.items():
            for i in range(0, len(binList)-1):
                for j in range(i+1, len(binList)):
                    # If these bins are not linked, link them
                    if not self.graph.has_edge(binList[i], binList[j]):
                        self.graph.add_edge(binList[i], binList[j], id={seqID})
                    
                    # Otherwise, modify the edge ID attribute
                    else:
                        self.graph[binList[i]][binList[j]]["id"].add(seqID)
    
    def prune(self, WEIGHT_CUTOFF=0.5):
        '''
        Curates the graph by removing edges that are below a certain weight threshold,
        removing sequence IDs that occur in many of the cut edges, and removes bins (nodes)
        that are left isolated after the pruning.
        
        Parameters:
            WEIGHT_CUTOFF -- a float greater than 0, and less than or equal to 1.0; this
                             value controls the weight threshold for which edges are cut,
                             wherein any edge with a weight less than (weightiest *
                             WEIGHT_CUTOFF) will be cut.
        '''
        binsToRemove = []
        for connectedBins in self.connected_components():
            connectedBins = list(connectedBins)
            
            # Skip over single bin components
            if len(connectedBins) == 1:
                continue
            
            # Obtain data for all edges in this component
            edgeData = []
            for edge in nx.edges(self.graph, connectedBins):
                edgeIDs = self[edge[0]][edge[1]]["id"]
                edgeData.append([
                    edge,
                    edgeIDs,
                    len(edgeIDs)
                ])
            
            # Cut any edges that are < (weightiest * WEIGHT_CUTOFF)
            weightiest = max([ ed[2] for ed in edgeData ])
            lowerBoundary = int(weightiest * WEIGHT_CUTOFF)
            
            cutEdges = {}
            for edge, edgeIDs, edgeWeight in edgeData:
                if edgeWeight < lowerBoundary:
                    self.graph.remove_edge(*edge)
                    for eID in edgeIDs:
                        cutEdges.setdefault(eID, 0)
                        cutEdges[eID] += 1
            
            # Remove the ID from all edges if relevant
            for seqID, numCuts in cutEdges.items():
                if numCuts >= lowerBoundary: # this cutoff should suffice as a scaling value
                    self.eliminations.add(seqID) # permanently eliminate from BINge clustering
                    
                    for edge, edgeIDs, edgeWeight in edgeData:
                        if seqID in edgeIDs and self.graph.has_edge(*edge):
                            self[edge[0]][edge[1]]["id"].remove(seqID)
            
            # Flag bins for removal after iteration
            for bin in connectedBins:
                if nx.is_isolate(self.graph, bin):
                    binsToRemove.append(bin)
        
        # Remove bins flagged for removal
        for bin in binsToRemove:
            self.graph.remove_node(bin)
    
    def connected_components(self):
        return nx.connected_components(self.graph)
    
    def __getitem__(self, key):
        return self.graph.__getitem__(key)
    
    def __len__(self):
        return len(self.graph)
    
    def __repr__(self):
        return "<BinGraph object;num_bins={0}>".format(
            len(self.graph)
        )

class BinBundle:
    '''
    Encapsulates a simple list of Bin objects for storage and manipulation when the searching
    of an IntervalTree (as provided by BinCollection) is not needed.
    
    The naming is arbitrary. A bundle is less organised than a collection I suppose?
    '''
    def __init__(self):
        self.bins = []
        self.eliminations = None
    
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
    
    def cluster_by_cooccurrence(self, VOTE_THRESHOLD = 0.66):
        '''
        This algorithm clusters sequences based on their co-occurrence in exon bins.
        Sequences which co-occur in bins at a certain frequency of their total number
        of bins that occur in are linked together. This linking occurs using a
        graph-based approach to modelling the relationships between the bins.
        
        Unlike link_by_sharing, this algorithm is guaranteed to produce a stable output
        and hence convergence is not a concern.
        
        Parameters:
            VOTE_THRESHOLD -- a float greater than 0, and less than or equal to 1.0.
                              This represents a ratio for which sequences must co-occur.
        Returns:
            clusterDict -- a dictionary where keys are integers from 0 -> n, and values
                           are sets of sequence IDs.
        '''
        assert 0 < VOTE_THRESHOLD <= 1.0, \
            "VOTE_THRESHOLD must be a value greater than 0, and less than or equal to 1"
        
        # Figure out which bins each sequence occur in
        occurrenceDict = {}
        for binIndex, bin in enumerate(self.bins):
            for seqID in bin.ids:
                occurrenceDict.setdefault(seqID, set())
                occurrenceDict[seqID].add(binIndex)
        
        # Initialise graph data structure to link sequences
        binGraph = nx.Graph()
        binGraph.add_nodes_from(occurrenceDict.keys())
        
        # Iterate through each sequence and link where appropriate
        for seqID, binIndices in occurrenceDict.items():
            # Count how many times other sequences occur in the same bin as this
            seqLinkDict = {}
            for binIndex in binIndices:
                bin = self.bins[binIndex]
                for otherSeqID in bin.ids:
                    if otherSeqID == seqID:
                        continue
                    
                    seqLinkDict.setdefault(otherSeqID, 0)
                    seqLinkDict[otherSeqID] += 1
            
            # See if any sequences meet the threshold for linking
            for otherSeqID, numOccurrence in seqLinkDict.items():
                if (numOccurrence / len(binIndices)) >= VOTE_THRESHOLD:
                    binGraph.add_edge(seqID, otherSeqID)
        
        # Create sequence bins based on the graph's connected components
        clusterDict = {}
        ongoingCount = 0
        for connectedSeqIDs in nx.connected_components(binGraph):
            clusterDict[ongoingCount] = set(connectedSeqIDs)
            ongoingCount += 1
        
        return clusterDict
    
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
    def create_from_graph(binGraph):
        '''
        Initializes a BinBundle object from a BinGraph object.
        
        Parameters:
            binGraph -- a BinGraph object.
        Returns:
            binBundle -- a BinBundle object.
        '''
        newBundle = BinBundle()
        for bin in binGraph.graph:
            newBundle.add(bin)
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
    
    @staticmethod
    def create_from_multiple_graphs(binGraphList):
        '''
        Initializes a BinBundle object from one or more BinGraph objects.
        The input value is expected to be a list of BinGraph objects.
        The result is a single BinBundle object with all the bins squished
        down into it.
        
        Parameters:
            binGraphList -- a list of BinGraph objects.
        Returns:
            binBundle -- a BinBundle object.
        '''
        assert isinstance(binGraphList, list), \
            "binGraphList must be a list of BinGraph objects"
        
        # Get all the eliminations sequence IDs
        allEliminations = set()
        for binGraph in binGraphList:
            allEliminations = allEliminations.union(binGraph.eliminations)
        
        # Generate the new bundle, removing any eliminated sequences
        newBundle = BinBundle()
        for binGraph in binGraphList:
            for bin in binGraph.graph:
                # Remove eliminated sequences from any bins
                toRemove = bin.ids.intersection(allEliminations)
                if len(toRemove) > 0:
                    bin.ids = bin.ids.difference(toRemove)
                
                # Add the bin to the bundle if it still contains IDs
                if len(bin.ids) > 0:
                    newBundle.add(bin)
        
        newBundle.eliminations = allEliminations
        return newBundle
    
    def __len__(self):
        return len(self.bins)
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinBundle object;num_bins={0}>".format(
            len(self)
        )
