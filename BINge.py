#! python3
# BINge.py
# BIN Genes for Expression analyses

# This script/program aims to provide an alternative to Corset
# for situations where a reference genome is available. However,
# the reference genome might be for a different subspecies than
# the one you're working with. Hence, it's not good enough to
# just map against the reference genome. You should map against a
# de novo transcriptome, but leverage the genomic information
# to group transcripts into genes. That's what this does.

import os, argparse, sys, queue, time, pickle # REMOVE AFTER TESTING
import networkx as nx

from intervaltree import IntervalTree, Interval
from pyfaidx import Fasta
from threading import Thread

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts import ZS_GFF3IO, ZS_ClustIO, ZS_SeqIO

# Define classes
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
        self.ids = self.ids.union(otherBin.ids)
    
    def __repr__(self):
        return (f"<Bin object;contig='{self.contig}';start={self.start};" +
                f"end={self.end};num_ids={len(self.ids)}"
        )

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

class GmapBinThread(Thread):
    '''
    This provides a modified Thread which allows for the output of bin_by_gmap
    to be stored locally within the object. Multi-threading this process is viable
    since each GMAP file should correspond to an entirely separate genome and hence
    there should be no overlap.
    '''
    def __init__(self, gmapFile, binCollection):
        Thread.__init__(self)
        
        self.gmapFile = gmapFile
        self.binCollection = binCollection
        self.novelBinCollection = None
        self.chimeras = None
    
    def run(self):
        self.novelBinCollection, self.chimeras = bin_by_gmap(self.gmapFile, self.binCollection)

class OutputWorkerThread(Thread):
    '''
    This provides a modified Thread which allows for the output worker to return a value
    indicating how many clusters it wrote to file. This will save a little bit of time later
    when we want to append more cluster results to the same file
    '''
    def __init__(self, outputQueue, outputFileName, clusterType):
        Thread.__init__(self)
        
        self.outputQueue = outputQueue
        self.outputFileName = outputFileName
        self.clusterType = clusterType
        self.numClusters = None
    
    def run(self):
        self.numClusters = output_worker(self.outputQueue, self.outputFileName, self.clusterType)

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.transcriptomeFile):
        print(f'I am unable to locate the transcriptome FASTA file ({args.transcriptomeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for annotationFile in args.annotationFiles:
        if not os.path.isfile(annotationFile):
            print(f'I am unable to locate the genome GFF3 file ({annotationFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for gmapFile in args.gmapFiles:
        if not os.path.isfile(gmapFile):
            print(f'I am unable to locate the GMAP GFF3 file ({gmapFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate the input files are logically sound
    if len(args.annotationFiles) != len(args.gmapFiles):
        print("Your genome annotation and GMAP files are incompatible!")
        print(f"I'm seeing {len(args.annotationFiles)} annotation files and {len(args.gmapFiles)} GMAP files")
        print("These numbers should be the same. You need to fix this up and try again.")
        quit()
    # Validate numeric inputs
    if args.threads < 1:
        print("--threads should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if args.convergenceIters < 1:
        print("--convergence_iters should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if not 0.0 < args.cdhitIdentity <= 1.0:
        print("--cdhit_identity should be given a value greater than zero, and equal or less than 1")
        print("Fix this and try again.")
        quit()
    if not 0.0 <= args.cdhitShortCov <= 1.0:
        print("--cdhit_shortcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not 0.0 <= args.cdhitLongCov <= 1.0:
        print("--cdhit_longcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not args.cdhitMem >= 100:
        print("--cdhit_mem should be given a value at least greater than 100 (megabytes)")
        print("Fix this and try again.")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def generate_bin_collections(annotationFiles):
    '''
    Receives a list of genome annotations in GFF3 format and parses these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        annotationFiles -- a list of strings pointing to GFF3 genome annotation files.
    Returns:
        collectionList -- a list containing BinCollections in the same order as the
                          input annotation files.
    '''
    # Parse each GFF3 into a bin collection structure
    collectionList = [] # by using a list, we keep genome bins separate to multi-thread later
    
    for annotFile in annotationFiles:
        binCollection = BinCollection()
        
        # Parse first as a GFF3 object
        gff3Obj = ZS_GFF3IO.GFF3(annotFile, strict_parse=True)
        for geneFeature in gff3Obj.types["gene"]:
            
            # Create a bin for this feature
            featureBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
            featureBin.add(geneFeature.ID)
            
            # See if this overlaps an existing bin
            binOverlap = binCollection.find(geneFeature.contig, geneFeature.start, geneFeature.end)
            
            # If not, add the new bin
            if len(binOverlap) == 0:
                binCollection.add(featureBin)
            
            # Otherwise...
            else:
                # ... merge the bins together
                for overlappingBin in binOverlap:
                    featureBin.merge(overlappingBin)
                
                # ... delete the overlapping bins
                for overlappingBin in binOverlap:
                    binCollection.delete(overlappingBin)
                
                # ... and add the new bin
                binCollection.add(featureBin)
        
        # Store the bin collection in our list for multi-threading later
        collectionList.append(binCollection)
    
    return collectionList

def _create_feature_from_sl(sl):
    # Extract details from sl
    contig, source, featureType, start, end, \
        score, strand, frame, attributes \
        = sl
    
    # Build feature with main details
    feature = ZS_GFF3IO.Feature()
    feature.add_attributes({
        "contig": contig, "source": source, "type": featureType,
        "start": int(start), "end": int(end), "coords": [int(start), int(end)],
        "score": score, "strand": strand, "frame": frame
    })
    
    # Add attributes
    splitAttributes = []
    for a in attributes.split("="):
        if ";" in a:
            splitAttributes += a.rsplit(";", maxsplit=1)
        else:
            splitAttributes.append(a)
    attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
    feature.add_attributes(attributesDict)
    
    return feature

def iterate_through_gff3(gff3File):
    '''
    Provides a simple iterator for a GFF3 file which yields GFF3
    features. Only works if the GFF3 file is sorted.
    
    Parameters:
        gff3File -- a string indicating the location of a GFF3 formatted
                    file that is sorted.
    '''
    thisFeature = []
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\t\r\n ").split("\t")
            
            # Skip comment / irrelevant lines
            if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                continue
            
            # Handle content lines
            else:
                featureType = sl[2]
                
                # Yield a completed gene feature
                if featureType == "gene" and thisFeature != []:
                    # Create the base gene and mRNA features
                    feature = _create_feature_from_sl(thisFeature[0])
                    mrnaFeature = _create_feature_from_sl(thisFeature[1])
                    feature.add_child(mrnaFeature)
                    
                    # Add subfeatures to base
                    for _sl in thisFeature[2:]:
                        _feature = _create_feature_from_sl(_sl)
                        mrnaFeature.add_child(_feature)
                    
                    # Reset our feature storage, and yield result now
                    thisFeature = []
                    yield feature
                
                # Build an ongoing feature
                thisFeature.append(sl)

def calculate_overlap_percentages(feature, bin):
    '''
    Feature should come from the GMAP parsing.
    Bin should come from a BinCollection.
    '''
    overlapLength = min(feature.end, bin.end) - max(feature.start, bin.start)
    featureOvlPct = overlapLength / (feature.end - feature.start + 1)
    binOvlPct = overlapLength / (bin.end - bin.start + 1)
    
    return featureOvlPct, binOvlPct

def novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT):
    '''
    Helper function pulled out to prevent code repetition. It will
    take a feature that has been determined to be "novel" and handle
    whether it should be binned with another novel bin or become
    a new bin itself.
    '''
    baseID = feature.ID.rsplit(".", maxsplit=1)[0]
    
    thisBin = Bin(feature.contig, feature.start, feature.end) # create a novel bin
    thisBin.add(baseID) # we don't want the .path# suffix attached to this bin
    
    # ... see if this overlaps a novel bin that it should join
    novelBinOverlap = novelBinCollection.find(feature.contig, feature.start, feature.end)
    binsToJoin = []
    for novelBin in novelBinOverlap:
        featureOvlPct, binOvlPct = calculate_overlap_percentages(feature, novelBin)
        shouldJoin = featureOvlPct >= NOVEL_BIN_OVL_PCT or binOvlPct >= NOVEL_BIN_OVL_PCT
        if shouldJoin:
            binsToJoin.append(novelBin)
    
    # Merge any bins that should join based on this novel sequence
    for novelBin in binsToJoin:
        thisBin.merge(novelBin)
        novelBinCollection.delete(novelBin)
    
    # Store this bin in the collection
    novelBinCollection.add(thisBin)

def bin_by_gmap(gmapFile, binCollection):
    '''
    Parameters:
        gmapFile -- a string pointing to the location of a GMAP GFF3 file
                    containing transcript alignments.
        binCollection -- a BinCollection object containing gene feature bins
                         generated from an official genome annotation.
    '''
    # Behavioural parameters (static for now, may change later)
    OKAY_COVERAGE = 94
    OKAY_IDENTITY = 92
    ##
    ALLOWED_COV_DIFF = 1
    ALLOWED_IDENT_DIFF = 0.1
    ##
    GENE_BIN_OVL_PCT = 0.5
    NOVEL_BIN_OVL_PCT = 0.5
    
    # Parse GMAP file and begin binning
    novelBinCollection = BinCollection() # holds onto novel bins we create during this process
    chimeras = set() # holds onto sequences we should filter if we didn't find a good path
    pathDict = {} # holds onto sequences we've already checked a path for
        
    for feature in iterate_through_gff3(gmapFile):
        mrnaFeature = feature.mRNA[0]
        
        # Get alignment statistics
        coverage, identity = float(mrnaFeature.coverage), float(mrnaFeature.identity)
        
        # Skip processing if the alignment sucks
        isGoodAlignment = coverage >= OKAY_COVERAGE and identity >= OKAY_IDENTITY
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
        
        # Now that we've checked if this path is good, we'll store it
        """This means on the next loop if there's another path for this gene, but it's worse
        than this current one, we'll just skip it"""
        pathDict.setdefault(baseID, [coverage, identity])
        
        # See if this overlaps an existing bin
        binOverlap = binCollection.find(feature.contig, feature.start, feature.end)
        
        # If it does not overlap an existing bin ...
        if len(binOverlap) == 0:
            # ... either add it to an existing novel bin or create a new novel bin
            novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT)
        
        # Exclude any transcripts that overlap multiple genes
        elif len(binOverlap) > 1:
            "We expect these occurrences to be chimeric transcripts which should be filtered"
            if baseID not in pathDict: # we don't want to filter a transcript if it's already
                chimeras.add(baseID)   # had a good path found for it
        
        # Compare to the existing bin
        else:
            geneBin = binOverlap[0]
            featureOvlPct, binOvlPct = calculate_overlap_percentages(feature, geneBin)
            shouldJoin = featureOvlPct >= GENE_BIN_OVL_PCT or binOvlPct >= GENE_BIN_OVL_PCT
            
            # If this should join the gene bin, do so now
            if shouldJoin:
                geneBin.add(baseID)
            
            # Otherwise, add it to an existing novel bin or create a new novel bin
            else:
                novel_binner(feature, novelBinCollection, NOVEL_BIN_OVL_PCT)
    
    return novelBinCollection, chimeras

def cluster_bin(bin, transcriptRecords):
    '''
    Parameters:
        bin -- a Bin object which should be clustered with CD-HIT.
        transcriptRecords -- a pyfaidx Fasta object containing all the sequences
                             we'll want to cluster from our bins.
    Returns:
        clusterDict -- a dictionary result of clustering with structure like:
                       {
                           0: [seqid1, seqid2, ...],
                           1: [ ... ],
                           ...
                       }
    '''
    # Bypass clustering if this bin has only one sequence
    if len(bin.ids) == 1:
        clusterDict = {0: list(bin.ids)}
    else:
        # Create a FASTA object for the sequences in this bin
        FASTA_obj = ZS_SeqIO.FASTA(None)
        for id in bin.ids:
            FastASeq_obj = ZS_SeqIO.FastASeq(id, str(transcriptRecords[id]))
            FASTA_obj.insert(0, FastASeq_obj)
        
        # Cluster it with CD-HIT at lax settings
        clusterer = ZS_ClustIO.CDHIT(FASTA_obj, "nucleotide")
        clusterer.identity = 0.8
        clusterer.set_shorter_cov_pct(0.2)
        clusterer.set_longer_cov_pct(0.0)
        clusterer.set_local()
        clusterer.get_cdhit_results(returnFASTA=False, returnClusters=True) # throw away the temporary file return value
        
        clusterDict = clusterer.resultClusters
    
    return clusterDict

def bin_clustering_worker(binQueue, outputQueue, transcriptRecords):
    '''
    Parameters:
        binQueue -- a queue.Queue object containing Bin objects for clustering.
        outputQueue -- a queue.Queue object that will receive the processed outputs
                       of this worker thread.
        transcriptRecords -- a pyfaidx Fasta object containing all the sequences
                             we'll want to cluster from our bins.
    '''
    while True:
        # Continue condition
        if binQueue.empty():
            time.sleep(0.5)
            continue
        
        # Grabbing condition
        bin = binQueue.get()
        
        # Exit condition
        if bin == None:
            binQueue.task_done()
            break
        
        # Perform work
        clusterDict = cluster_bin(bin, transcriptRecords)
        
        # Store result in output queue
        outputQueue.put(clusterDict)
        
        # Mark work completion
        binQueue.task_done()

def output_worker(outputQueue, outputFileName, clusterType):
    '''
    Parameters:
        outputQueue -- a queue.Queue object that receives outputs from worker threads.
        outputFileName -- a string indicating the location to write results to.
        clusterType -- a string to put in an output column indicating what source of
                       evidence led to this cluster's formation; usually, this will
                       have the value of "binned" to indicate that it was binned via
                       GMAP alignment, or "unbinned" to indicate that it has been
                       binned solely via CD-HIT clustering.
    '''
    clusterNum = 1
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines as they become available
        while True:
            # Continue condition
            if outputQueue.empty():
                time.sleep(0.5)
                continue
            
            # Grabbing condition
            clusterDict = outputQueue.get()
            
            # Exit condition
            if clusterDict == None:
                outputQueue.task_done()
                break
            
            # Perform work
            for clusterIDs in clusterDict.values():
                for seqID in clusterIDs:
                    fileOut.write(f"{clusterNum}\t{seqID}\t{clusterType}\n")
                clusterNum += 1
            
            # Mark work completion
            outputQueue.task_done()
    return clusterNum - 1 # -1 to return the amount of clusters we actually wrote

def bin_self_linker(binCollection):
    '''
    Receives a BinCollection that has potentially been created through multiple
    genome, multiple subspecies transcriptomics. It will attempt to merge bins
    that are "equivalent" across the genomes using a graph-based approach to
    modelling the relationships between the bins. The result is a BinCollection
    that unifies all input data and is appropriate for sequence clustering.
    
    Note: this function is not guaranteed to produce a stable output i.e., if
    you feed its result back in, you may get different results. It should
    eventually converge on a stable solution, but there might theoretically
    be an edge case where that never happens.
    
    Parameters:
        binCollection -- a BinCollection object containing an assortment
                         of bins from different genomes and/or genes that
                         are (at a sequence level) indistinguishable and
                         hence confounding for DGE.
    Returns:
        linkedBinCollection -- a BinCollection object where bins have been
                               merged where deemed appropriate.
    '''
    VOTE_THRESHOLD = 0.5
    
    # Format bins into a dictionary
    "It's not ideal for memory use but I need instant lookup and consistent ordering"
    binDict = {}
    binCounter = 0
    for bin in binCollection:
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
            newBin.union(binDict[index].ids)
        linkedBinCollection.add(newBin)
    
    return linkedBinCollection

def iterative_bin_self_linking(binCollection, convergenceIters):
    '''
    Links a bin to itself iteratively until it converges or
    convergenceIters is reached.
    
    Parameters:
        binCollection -- a BinCollection object
        convergenceIters -- an integer providing a maximum limit for
                            the amount of iterations that may occur.
    '''
    prevCollectionCount = len(binCollection.bins)
    for _ in range(convergenceIters):
        binCollection = bin_self_linker(binCollection)
        if len(binCollection.bins) == prevCollectionCount:
            break
        prevCollectionCount = len(binCollection.bins)
    return binCollection

def multithread_bin_cluster(binCollectionList, threads, transcriptRecords,
                            outputFileName, clusterType="binned"):
    '''
    This code has been pulled out into a function solely to ensure the main() function
    is cleaner to read.
    
    What it does is set up a queue to store bins. Worker threads will grab bins to
    process with CD-HIT. They will put their predicted cluster dictionary into the output
    worker which will write it to file.
    
    We use this CD-HIT clustering of bins because each bin is only guaranteed to contain
    sequences that generally align well over a predicted genomic region. There's no
    guarantee that they are isoforms of the same gene, or even share any exonic content.
    This clustering step is very lax and designed just to separate out overlapping
    genes (which aren't actually the _same_ gene).
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        threads -- an integer indicating how many threads to run
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
        outputFileName -- a string indicating the location to write output clustering
                          results to.
        clusterType -- a string to put in an output column indicating what source of
                       evidence led to this cluster's formation; usually, this will
                       have the value of "binned" to indicate that it was binned via
                       GMAP alignment, or "unbinned" to indicate that it has been
                       binned solely via CD-HIT clustering.
    Returns:

    '''
    # Set up queueing system for multi-threading
    workerQueue = queue.Queue(maxsize=int(threads * 10)) # allow queue to scale with threads
    outputQueue = queue.Queue(maxsize=int(threads * 100))
    
    # Start up threads for clustering of bins
    for _ in range(threads):
        worker = Thread(
            target=bin_clustering_worker,
            args=(workerQueue, outputQueue, transcriptRecords))
        worker.setDaemon(True)
        worker.start()
    
    outputWorker = OutputWorkerThread(outputQueue, outputFileName, clusterType)
    outputWorker.setDaemon(True)
    outputWorker.start()
    
    # Put bins in queue for worker threads
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            workerQueue.put(bin)
    
    # Close up shop on the threading structures
    for i in range(threads):
        workerQueue.put(None) # this is a marker for the worker threads to stop
    workerQueue.join()
    
    outputQueue.put(None) # marker for the output thead to stop
    outputQueue.join()
    
    return outputWorker.numClusters

def get_unbinned_sequence_ids(binCollectionList, transcriptRecords):
    '''
    Compares one or more BinCollection objects against the transcript sequences
    to see if any sequences indicated in transcriptRecords do not exist in
    any Bins.
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    '''
    binnedIDs = []
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            binnedIDs.extend(bin.ids)
    binnedIDs = set(binnedIDs)
    
    unbinnedIDs = []
    for record in transcriptRecords:
        if record.name not in binnedIDs:
            unbinnedIDs.append(record.name)
    return unbinnedIDs

def cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, threads, mem,
                               IDENTITY=0.85, SHORTER_COV_PCT=0.6, LONGER_COV_PCT=0.3):
    '''
    Runs CD-HIT on the unbinned sequences in order to assign them to a cluster.
    It's not ideal, but the alternative is to exclude these sequences which may
    not be in line with the user's aims. By noting the source of clustering in
    the output file, they can decide if they'd like to exclude these or not.
    
    Note: the default CD-HIT parameter values have been derived from some limited
    testing which suggests that these values may work optimally to cluster sequences
    into "genes" or the closest thing we have to them without proper genomic evidence.
    
    Parameters:
        unbinnedIDs -- a list containg sequence IDs that exist in transcriptRecords
                       which should be clustered with CD-HIT.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences.
        threads -- an integer value indicating how many threads to run CD-HIT with.
        mem -- an integer value indicating how many megabytes of memory to run CD-HIT with.
        IDENTITY -- a float value indicating what identity value to run CD-HIT with.
        SHORTER_COV_PCT -- a float value setting the -aS parameter of CD-HIT.
        LONGER_COV_PCT -- a float value setting the -aL parameter of CD-HIT.
    '''
    # Generate a temporary FASTA file containing unbinned transcripts
    tmpFileName = "tmp_BINge_unbinned_{0}.fasta".format(
        ZS_SeqIO.Conversion.get_hash_for_input_sequences(str(transcriptRecords))
    )
    with open(tmpFileName, "w") as fileOut:
        for seqID in unbinnedIDs:
            record = transcriptRecords[seqID]
            fileOut.write(f">{record.name}\n{str(record)}\n")
    
    # Cluster the unbinned transcripts
    clusterer = ZS_ClustIO.CDHIT(tmpFileName, "nucleotide")
    clusterer.identity = IDENTITY
    clusterer.set_shorter_cov_pct(SHORTER_COV_PCT)
    clusterer.set_longer_cov_pct(LONGER_COV_PCT)
    clusterer.threads = threads
    clusterer.mem = mem
    clusterer.get_cdhit_results(returnFASTA=False, returnClusters=True)
    
    # Clean up temporary file(s)
    os.unlink(tmpFileName)
    
    # Return cluster dictionary results
    return clusterer.resultClusters

## Main
def main():
    # User input
    usage = """%(prog)s (BIN Genes for Expression analyses) is a program which bins
    de novo-assembled transcripts together on the basis of reference genome alignments.
    This might be necessary when working with multiple subspecies that are expected to
    diverge only slightly from the reference organism. By binning like this, each subspecies
    can have its gene counts compared fairly during DGE.
    
    Note: For each annotation GFF3 (-ga) you should provide a matching GMAP GFF3 alignment
    file (-gm). These values should be ordered equivalently.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="transcriptomeFile",
                   required=True,
                   help="Input transcriptome FASTA file (mRNA or CDS).")
    p.add_argument("-ga", dest="annotationFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more genome annotation (GFF3) files")
    p.add_argument("-gm", dest="gmapFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more GMAP (GFF3) files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    # Optional
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many threads to run when multithreading
                   is available (default==1)""",
                   default=1)
    p.add_argument("--convergence_iters", dest="convergenceIters",
                   required=False,
                   type=int,
                   help="""Optionally, specify a maximum number of iterations allowed for
                   bin convergence to be achieved (default==5); in most cases results will
                   converge in fewer than 5 iterations, so setting a maximum acts merely as
                   a safeguard against edge cases I have no reason to believe will ever
                   happen.""",
                   default=5)
    p.add_argument("--cdhit_identity", dest="cdhitIdentity",
                   required=False,
                   type=float,
                   help="""Optionally, when clustering unbinned sequences, specify what
                   identity value (-c) to provide CD-HIT (default==0.85)""",
                   default=0.85)
    p.add_argument("--cdhit_shortcov", dest="cdhitShortCov",
                   required=False,
                   type=float,
                   help="""Optionally, when clustering unbinned sequences, specify what
                   -aS parameter to provide CD-HIT (default==0.6)""",
                   default=0.6)
    p.add_argument("--cdhit_longcov", dest="cdhitLongCov",
                   required=False,
                   type=float,
                   help="""Optionally, when clustering unbinned sequences, specify what
                   -aL parameter to provide CD-HIT (default==0.3)""",
                   default=0.3)
    p.add_argument("--cdhit_mem", dest="cdhitMem",
                   required=False,
                   type=int,
                   help="""Optionally, when clustering unbinned sequences, specify how
                   many megabytes of memory to provide CD-HIT (default==5000)""",
                   default=5000)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each GFF3 into a bin collection structure
    collectionList = generate_bin_collections(args.annotationFiles) # keep genome bins separate to multi-thread later
    
    # Parse GMAP alignments into our bin collection with multiple threads
    novelBinCollection = BinCollection()
    chimeras = set()
    for i in range(0, len(args.gmapFiles), args.threads): # only process n (threads) files at a time
        processing = []
        for x in range(args.threads): # begin processing n files
            if i+x < len(args.gmapFiles): # parent loop may excess if n > the number of GMAP files
                gmapFile = args.gmapFiles[i+x]
                binCollection = collectionList[i+x]
                
                workerThread = GmapBinThread(gmapFile, binCollection)
                processing.append(workerThread)
                workerThread.start()
        
        # Gather results
        for workerThread in processing:
            # Wait for thread to end
            workerThread.join()
            
            # Grab the new outputs from this thread
            "Each thread modifies a BinCollection part of collectionList directly"
            threadNovelBinCollection, threadChimeras = workerThread.novelBinCollection, \
                workerThread.chimeras
            
            # Merge them
            novelBinCollection.merge(threadNovelBinCollection)
            chimeras = chimeras.union(threadChimeras)
    
    # Merge gene bins together
    binCollection = collectionList[0]
    for i in range(1, len(collectionList)):
        binCollection.merge(collectionList[i])
    
    # Link and merge bins across genomes / across gene copies
    """Usually linking will unify multiple genomes together, but it may detect
    bins of identical gene copies and link them together which is reasonable
    since these would confound DGE to keep separate anyway"""
    
    binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
    novelBinCollection = iterative_bin_self_linking(novelBinCollection, args.convergenceIters)
    
    # Merge bin collections together and re-link bins
    """A novel bin in one genome may just be because of an absence in the annotation.
    Such an absence may not exist in another genome, and hence it will have a gene bin.
    These should be merged together to prevent having redundant bins."""
    
    binCollection.merge(novelBinCollection)
    binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
    
    # Resolve chimeras
    ## For implementation only if chimeras are found in a dataset I can test on
    
    # Load transcripts into memory for quick access
    transcriptRecords = Fasta(args.transcriptomeFile)
    
    # Cluster each bin with CD-HIT to separate overlapping genes
    numClusters = multithread_bin_cluster([binCollection],
                                          args.threads, transcriptRecords,
                                          args.outputFileName, clusterType="binned")
    
    # Cluster remaining unbinned sequences
    unbinnedIDs = get_unbinned_sequence_ids([binCollection], transcriptRecords)
    unbinnedClusterDict = cluster_unbinned_sequences(
        unbinnedIDs, transcriptRecords, args.threads, args.cdhitMem,
        args.cdhitIdentity, args.cdhitShortCov, args.cdhitLongCov)
    
    ## TESTING ##
    with open("BINge_testing.pkl", "wb") as fileOut:
        pickle.dump([binCollection, numClusters, unbinnedClusterDict], fileOut)
    ## END TESTING ##
    
    # Write output of clustering to file
    with open(args.outputFileName, "a") as fileOut:
        for clusterNum, clusterIDs in unbinnedClusterDict.items():
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum+numClusters+1}\t{seqID}\tunbinned\n") # clusterType = "unbinned"
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
