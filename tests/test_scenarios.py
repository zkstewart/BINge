#! python3

import os, sys, unittest

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_GFF3IO import GFF3
from modules.bins import BinCollection, Bin
from modules.thread_workers import GmapBinThread
from modules.bin_handling import generate_bin_collections, populate_bin_collections, \
    multithread_bin_splitter, iterative_bin_self_linking

###

def _generate_bin_collections(gff3Files):
    '''
    A re-implement with the same logic just minus some of the setup.
    '''
    # Generate and seed bin collections
    collectionList = []
    
    for gff3File in gff3Files:
        binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if gff3File != None:
            gff3Obj = GFF3(gff3File, strict_parse=False)
            for geneFeature in gff3Obj.types["gene"]:
                
                # Create a bin for this feature
                featureBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
                try:
                    featureBin.add(geneFeature.ID, Bin.format_exons_from_gff3_feature(geneFeature))
                except:
                    "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                    continue
                
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

def _populate_bin_collections(collectionList, gmapFiles, threads=1, gmapIdentity=0.95):
    '''
    A re-implement with the same logic just minus some of the setup.
    '''
    # Establish data structures for holding onto results
    novelBinCollection = BinCollection() # consolidated after each thread
    
    novelCollections = [BinCollection() for _ in range(len(collectionList))] # consolidated within each thread
    multiOverlaps = [[] for _ in range(len(collectionList))] # consolidated within threads
    
    # Establish lists for feeding data into threads
    threadData = []
    for i in range(len(collectionList)):
        # Get all GMAP files associated with this genome
        thisGmapFiles = [gmapFiles[i]]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[i]
        thisNovelCollection = novelCollections[i]
        thisMultiOverlap = multiOverlaps[i]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection, thisNovelCollection, thisMultiOverlap])
    
    # Start up threads
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection, thisNovelCollection, \
                    thisMultiOverlap = threadData[i+x]
                
                populateWorkerThread = GmapBinThread(thisGmapFiles, thisBinCollection,
                                                     thisNovelCollection, thisMultiOverlap,
                                                     gmapIdentity)
                
                processing.append(populateWorkerThread)
                populateWorkerThread.start()
        
        # Gather results
        for populateWorkerThread in processing:
            # Wait for thread to end
            populateWorkerThread.join()
            
            # Grab the new outputs from this thread
            """Each thread modifies a BinCollection part of collectionList directly,
            as well as a list in multiOverlaps"""
            threadNovelBinCollection = populateWorkerThread.novelBinCollection
            
            # Merge them via flattening
            novelBinCollection.flatten(threadNovelBinCollection)
    
    return novelBinCollection, multiOverlaps

###

# Specify data locations
dataDir = os.path.join(os.getcwd(), "data")

# Define pipeline for running the core BINge operations
def binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5, gmapIdentity=None):
    if gmapIdentity == None:
        gmapIdentity = 0.95
    
    novelBinCollection, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                 threads, gmapIdentity)
    
    for i in range(len(binCollectionList)):
        binCollection = binCollectionList[i].fix_fragments(multiOverlaps[i])
        binCollectionList[i] = binCollection
    
    binCollection = binCollectionList[0]
    for i in range(1, len(binCollectionList)):
        binCollection.merge(binCollectionList[i])
    
    binCollection = multithread_bin_splitter(binCollection, threads)
    
    binCollection = iterative_bin_self_linking(binCollection, convergenceIters)
    novelBinCollection = iterative_bin_self_linking(novelBinCollection, convergenceIters)
    
    binCollection.merge(novelBinCollection)
    binCollection = iterative_bin_self_linking(binCollection, convergenceIters)
    
    return binCollection, novelBinCollection

# Define test helper functions
def get_binCollection_in_range(gff3File, contig, start, end):
    binCollectionList = _generate_bin_collections([gff3File])
    binCollection = binCollectionList[0]
    
    newBinCollection = BinCollection()
    for bin in binCollection.find(contig, start, end):
        newBinCollection.add(bin)
    return newBinCollection

# Define test scenarios
class TestNormal(unittest.TestCase):
    def test_normal_1(self):
        '''
        In this test, we expect 1 ID since the base name of the GMAP path will be
        identical to the pre-seeded ID from annotation.gff3. This is expected behaviour
        since it 1) prevents redundancy, and 2) GMAP's idea of the model coordinates might
        differ to the annotation; we trust the annotation more then GMAP in this scenario.
        '''
        # Arrange
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 1, 100)
        binCollectionList = [binCollection]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollection, novelBinCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
        binOverlap = binCollection.find("contig1", 1, 50)
        binOverlap2 = binCollection.find("contig1", 51, 100)
        
        # Assert
        self.assertEqual(len(binOverlap), 1, "Should contain 1 bin")
        self.assertEqual(len(binOverlap2), 1, "Should contain 1 bin")
        
        self.assertEqual(len(binOverlap[0].ids), 1, "Should contain 1 ID")
        self.assertEqual(len(binOverlap2[0].ids), 1, "Should contain 1 ID")
    
    def test_normal_2(self):
        '''
        In this test, we change the annotation bin IDs to allow the GMAP paths to be added
        to them.
        '''
        # Arrange
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 1, 100)
        for i, bin in enumerate(binCollection):
            bin.data.ids = {f"example{i+1}"}
            if i == 0:
                bin.data.exons = {f"example{i+1}":[[1, 50]]}
            else:
                bin.data.exons = {f"example{i+1}":[[51, 100]]}
        
        binCollectionList = [binCollection]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollection, novelBinCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
        binOverlap = binCollection.find("contig1", 1, 50)
        binOverlap2 = binCollection.find("contig1", 51, 100)
        
        # Assert
        self.assertEqual(len(binOverlap), 1, "Should contain 1 bin")
        self.assertEqual(len(binOverlap2), 1, "Should contain 1 bin")
        
        self.assertEqual(len(binOverlap[0].ids), 2, "Should contain 2 IDs")
        self.assertEqual(len(binOverlap2[0].ids), 2, "Should contain 2 IDs")

class TestBinSplitter(unittest.TestCase):
    def test_splitter_normal(self):
        '''
        This test should result in no splitting occurring
        '''
        # Arrange
        threads = 1
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 1, 100)
        for i, bin in enumerate(binCollection):
            bin.data.ids = {f"example{i+1}"}
            if i == 0:
                bin.data.exons = {f"example{i+1}": [[1, 50]]}
            else:
                bin.data.exons = {f"example{i+1}": [[51, 100]]}
        
        binCollectionList = [binCollection]
        gmapFiles = [os.path.join(dataDir, "gmap.gff3")]
        
        novelBinCollection, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads, 0.95)
        
        binCollection = binCollectionList[0]
        for i in range(1, len(binCollectionList)):
            binCollection.merge(binCollectionList[i])
        
        # Act
        binCollection = multithread_bin_splitter(binCollection, threads)
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should contain 2 bins")
    
    def test_splitter_nested(self):
        '''
        This test should result in splitting of the nested features from the
        parent gene feature.
        '''
        # Arrange
        threads = 1
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 210, 300)
        origNumBins = len(binCollection)
        
        # Act
        binCollection = multithread_bin_splitter(binCollection, threads)
        newNumBins = len(binCollection)
        
        # Assert
        self.assertEqual(origNumBins, 1, "Should contain 1 bin")
        self.assertEqual(newNumBins, 2, "Should contain 2 bins")
    
    def test_splitter_overlapping(self):
        '''
        This test should result in splitting of the overlapping features since their
        length of overlap is minimal
        '''
        # Arrange
        threads = 1
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 310, 400)
        origNumBins = len(binCollection)
        origNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        # Act
        binCollection = multithread_bin_splitter(binCollection, threads)
        newNumBins = len(binCollection)
        newNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        # Assert
        self.assertEqual(origNumBins, 1, "Should contain 1 bin")
        self.assertEqual(newNumBins, 2, "Should contain 2 bins")
        
        self.assertEqual(origNumIDs, [4], "Should contain 1 bin with 4 IDs")
        self.assertEqual(newNumIDs, [2, 2], "Should contain 2 bins with 2 IDs each")

class TestFragmentMerger(unittest.TestCase):
    def test_fragment_merger(self):
        '''
        This test should result in the fragmented gene bins being merged together
        '''
        # Arrange
        gmapFile = os.path.join(dataDir, "annotation_fragments.gff3")
        binCollectionList = _generate_bin_collections([gmapFile])
        binCollection = binCollectionList[0]
        gmapFiles = [os.path.join(dataDir, "gmap_fragments.gff3")]
        
        # Act
        origNumBins = len(binCollection)
        origNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        binCollection, novelBinCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
        
        newNumBins = len(binCollection)
        newNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        # Assert
        self.assertEqual(origNumBins, 2, "Should contain 2 bins")
        self.assertEqual(newNumBins, 1, "Should contain 1 bin")
        
        self.assertEqual(origNumIDs, [1, 1], "Should contain 2 bins with 1 ID each")
        self.assertEqual(newNumIDs, [4], "Should contain 1 bin with 4 IDs")

if __name__ == '__main__':
    unittest.main()

threads=1
convergenceIters=5