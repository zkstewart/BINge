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

import os, argparse, sys, queue, time
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
        #self.links = [] # list of other Bin objects this one connects to
    
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
        self.ids.add(idValue)
    
    def union(self, idList):
        self.ids = self.ids.union(idList)
    
    def add_link(self, binValue):
        self.links.append(binValue)
    
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
    
    def __iter__(self):
        return iter(self.bins)
    
    def __repr__(self):
        return "<BinCollection object;num_bins={0}>".format(
            len(self.bins)
        )

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
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

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
    thisBin = Bin(feature.contig, feature.start, feature.end) # create a novel bin
    thisBin.add(feature.ID)
    
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
    ALLOWED_COV_DIFF = 2
    ALLOWED_IDENT_DIFF = 0.5
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
        clusterer.get_cdhit_results(returnClusters=True) # throw away the temporary file return value
        
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

def output_worker(outputQueue, outputFileName):
    '''
    Parameters:
        outputQueue -- a queue.Queue object that receives outputs from worker threads.
        outputFileName -- a string indicating the location to write results to.
    '''
    clusterNum = 1
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\n")
        
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
                    fileOut.write(f"{clusterNum}\t{seqID}\n")
                clusterNum += 1
            
            # Mark work completion
            outputQueue.task_done()

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
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each GFF3 into a bin collection structure
    binCollection = BinCollection()
    for annotFile in args.annotationFiles:
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
    gff3Obj = None # help garbage collection and reduce memory footprint (I think?)
    
    # For each GMAP file, add sequences to bins
    for gmapFile in args.gmapFiles:
        novelBinCollection, chimeras = bin_by_gmap(gmapFile, binCollection)
    
    # Resolve chimeras
    ## For implementation only if chimeras are found in dataset
    
    # Load transcripts into memory for quick access
    transcriptRecords = Fasta(args.transcriptomeFile)
    
    # Set up queueing system for multi-threading
    workerQueue = queue.Queue(maxsize=50)
    outputQueue = queue.Queue(maxsize=1000)
    
    # Start up threads for clustering of bins
    for _ in range(args.threads):
        worker = Thread(
            target=bin_clustering_worker,
            args=(workerQueue, outputQueue, transcriptRecords))
        worker.setDaemon(True)
        worker.start()
    
    outputWorker = Thread(target=output_worker,
                          args=(outputQueue, args.outputFileName))
    outputWorker.setDaemon(True)
    outputWorker.start()
    
    # Put bins in queue for worker threads
    for collection in [binCollection, novelBinCollection]:
        for interval in collection:
            bin = interval.data
            workerQueue.put(bin)
    
    # Close up shop on the threading structures
    for i in range(args.threads):
        workerQueue.put(None) # this is a marker for the worker threads to stop
    workerQueue.join()
    
    outputQueue.put(None) # marker for the output thead to stop
    outputQueue.join()
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
