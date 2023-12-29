#! python3

import os, sys, unittest, time
import networkx as nx
from multiprocessing import Pipe

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_GFF3IO import GFF3
from modules.bins import BinCollection, Bin
from modules.gff3_handling import iterate_through_gff3
from modules.bin_handling import generate_bin_collections, populate_bin_collections, \
    multithread_bin_splitter, iterative_bin_self_linking
from test_scenarios import get_binCollection_in_range, binge_runner
from modules.thread_workers import GmapBinThread

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
                    for mrnaFeature in geneFeature.mRNA:
                        featureBin.add(mrnaFeature.ID, Bin.format_exons_from_gff3_feature(mrnaFeature))
                    #featureBin.add(geneFeature.ID, Bin.format_exons_from_gff3_feature(geneFeature))
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
    multiOverlaps = [[] for _ in range(len(collectionList))] # consolidated within threads
    
    # Establish lists for feeding data into threads
    threadData = []
    for i in range(len(collectionList)):
        # Get all GMAP files associated with this genome
        thisGmapFiles = [gmapFiles[i]]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[i]
        thisMultiOverlap = multiOverlaps[i]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection, thisMultiOverlap])
    
    # Start up threads
    receivers = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection, \
                    thisMultiOverlap = threadData[i+x]
                thisReceiver, thisSender = Pipe()
                
                populateWorkerThread = GmapBinThread(thisGmapFiles, thisBinCollection,
                                                     thisMultiOverlap, thisSender,
                                                     gmapIdentity)
                
                processing.append(populateWorkerThread)
                populateWorkerThread.start()
                receivers.append(thisReceiver)
        
        # Gather results
        for populateWorkerThread in processing:
            # Wait for thread to end
            populateWorkerThread.join()
    
    # Gather results
    resultBinCollection, resultMultiOverlaps = [], []
    for receiver in receivers:
        binCollection, multiOverlap = receiver.recv()
        resultBinCollection.append(binCollection)
        resultMultiOverlaps.append(multiOverlap)
    
    return resultBinCollection, resultMultiOverlaps

###

# Specify data locations
dataDir = os.path.join(os.getcwd(), "data")

# Define unit tests
class TestBinCollection(unittest.TestCase):
    def test_binCollection_init(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "gmap_normal.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0]
        
        # Act
        ## N/A
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should be 2")
    
    def test_binCollection_find(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0]
        
        # Act
        binOverlap = binCollection.find("contig1", 1, 50)
        binOverlap2 = binCollection.find("contig1", 1, 100)
        
        # Assert
        self.assertEqual(len(binOverlap), 1, "Should be 1")
        
        self.assertEqual(binOverlap[0].start, 1, "Should be 1")
        self.assertEqual(binOverlap[0].end, 50, "Should be 50")
        self.assertEqual(binOverlap[0].contig, "contig1", "Should be contig1")
        
        self.assertEqual(len(binOverlap[0].ids), 1, "Should contain 1")
        self.assertIn("contig1.1.mrna", binOverlap[0].ids, "Should contain contig1.1.mrna")
        
        self.assertEqual(len(binOverlap[0].exons), 1, "Should contain 1")
        self.assertIn("contig1.1.mrna", binOverlap[0].exons, "Should contain contig1.1.mrna")
        self.assertEqual(binOverlap[0].exons["contig1.1.mrna"], [[1, 50]], "Should be [[1, 50]]")
        
        self.assertEqual(len(binOverlap2), 2, "Should be 2")
    
    def test_binCollection_delete(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList1 = _generate_bin_collections(gff3Files)
        binCollectionList2 = _generate_bin_collections(gff3Files)
        
        binCollection1 = binCollectionList1[0]
        binCollection2 = binCollectionList2[0]
        
        binOverlap1 = binCollection1.find("contig1", 1, 50)
        binOverlap2 = binCollection2.find("contig1", 51, 100)
        
        # Act
        beforeLength1 = len(binCollection1)
        beforeLength2 = len(binCollection2)
        
        binCollection1.delete(binOverlap1[0])
        binCollection2.delete(binOverlap2[0])
        
        afterLength1 = len(binCollection1)
        afterLength2 = len(binCollection2)
        
        # Assert
        self.assertEqual(beforeLength1, afterLength1 + 1, f"Should be {beforeLength1}")
        self.assertEqual(beforeLength2, afterLength2 + 1, f"Should be {beforeLength1}")
        
        for bin in binCollection1:
            self.assertNotEqual(bin.data.start, 1, f"Should not have start of 1")
        for bin in binCollection2:
            self.assertNotEqual(bin.data.start, 51, f"Should not have start of 51")

class TestBin(unittest.TestCase):
    def test_bin_add(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("contig1", 1, 100)
        bin1 = binOverlap[0]
        bin2 = binOverlap[1]
        
        newBinCollection = BinCollection()
        
        # Act 1
        newBinCollection.add(bin1)
        
        # Assert 1
        self.assertEqual(len(newBinCollection), 1, "Should contain 1 bin")
        
        # Act 2
        newBinCollection.add(bin2)
        
        # Assert 2
        self.assertEqual(len(newBinCollection), 2, "Should contain 2 bins")
    
    def test_bin_selfmerge(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("contig1", 1, 100)
        bin1 = binOverlap[0]
        bin2 = binOverlap[1]
        
        # Act
        bin1.merge(bin1)
        bin2.merge(bin2)
        
        # Assert
        self.assertEqual(bin2.start, 51, "Should be 51")
        self.assertEqual(bin2.end, 100, "Should be 100")
        self.assertEqual(bin2.contig, "contig1", "Should be contig1")
        
        self.assertEqual(len(bin1.ids), 1, "Should contain 1 id")
        self.assertEqual(len(bin2.ids), 1, "Should contain 1 id")
        
        self.assertIn("contig1.2.mrna", bin2.ids, "Should contain contig1.2.mrna")
    
    def test_bin_merge(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("contig1", 1, 100)
        bin1 = binOverlap[0]
        bin2 = binOverlap[1]
        
        # Act
        bin1.merge(bin2)
        
        # Assert
        self.assertEqual(bin1.start, 1, "Should be 1")
        self.assertEqual(bin1.end, 100, "Should be 100")
        self.assertEqual(bin1.contig, "contig1", "Should be contig1")
        
        self.assertEqual(len(bin1.ids), 2, "Should contain 2 ids")
        
        self.assertIn("contig1.1.mrna", bin1.ids, "Should contain contig1.1.mrna")   
        self.assertIn("contig1.2.mrna", bin1.ids, "Should contain contig1.2.mrna")

class TestGff3Iterate(unittest.TestCase):
    def test_iterate_through_empty_file(self):
        # Arrange
        gmapFile = os.path.join(dataDir, "empty_file")
        
        # Act and assert
        try:
            for feature in iterate_through_gff3(gmapFile):
                pass
        except:
            self.assertTrue(True, "Should error out here")
    
    def test_iterate_through_gff3(self):
        # Arrange
        gmapFile = os.path.join(dataDir, "gmap_normal.gff3")
        
        # Act
        numFeatures = 0
        mrnaIDs = []
        for feature in iterate_through_gff3(gmapFile):
            numFeatures += 1
            mrnaFeature = feature.mRNA[0]
            mrnaIDs.append(mrnaFeature.ID)
        
        # Assert
        self.assertEqual(numFeatures, 2, "Should contain 2 bins")
        self.assertEqual(mrnaIDs, ['contig1.1.mrna1', 'contig1.2.mrna1'],
                         "Should be ['contig1.1.mrna1', 'contig1.2.mrna1']")

class TestNovelPopulate(unittest.TestCase):
    def test_populate_reject_novel(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                  threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 0, "Should contain 0 bins")
    
    def test_populate_annotation_only(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 6, "Should contain 6 bins")
    
    def test_populate_novel_only(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should contain 2 bins")
    
    def test_populate_novel_merge(self):
        "I think this test isn't really testing what it's meant to anymore"
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 6, "Should contain 6 bins")
        for bin in binCollection:
            if bin.data.start == 1 or bin.data.start == 51:
                self.assertEqual(len(bin.data.ids), 2, f"Should have 2 IDs in bin")

class TestThreadExceptions(unittest.TestCase):
    def test_BinSplitWorker_has_exception(self):
        # Arrange
        binCollection = get_binCollection_in_range(os.path.join(dataDir, "annotation.gff3"),
                                                   "contig1", 1, 100)
        for i, bin in enumerate(binCollection):
            bin.data.ids = {f"example{i+1}"}
        
        binCollectionList = [binCollection]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        try:
            binCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
            self.assertTrue(False, "Not having exception is a test fail")
        # Assert
        except KeyError as e:
            self.assertEqual(eval(str(e)), "example1", "Should have exception at example1")
    
    def test_BinSplitWorker_not_has_exception(self):
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
        binCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
        
        # Assert
        self.assertTrue(True, "Not having exception is a test pass")

class TestNetworkX(unittest.TestCase):
    def test_connected_component_numbers(self):
        # Arrange
        binCollection = BinCollection()
        bin1 = Bin("contig1", 1, 100)
        bin1.ids = {"contig1.1"}
        bin1.exons = {"contig1.1": [[1, 100]]}
        
        bin2 = Bin("contig1", 99, 200)
        bin2.ids = {"contig1.2"}
        bin2.exons = {"contig1.2": [[99, 100200]]}
        
        bin3 = Bin("contig1", 199, 300)
        bin3.ids = {"contig1.3"}
        bin3.exons = {"contig1.3": [[199, 300]]}
        
        binCollection.add(bin1)
        binCollection.add(bin2)
        binCollection.add(bin3)
        
        # Act
        binDict = {}
        for bin in binCollection:
            binDict[hash(bin.data)] = bin.data

        binGraph = nx.Graph()
        binGraph.add_nodes_from(binDict.keys())
        
        binHashes = list(binDict.keys())
        for i in range(0, len(binHashes)-1):
            for x in range(i, len(binHashes)):
                binGraph.add_edge(binHashes[i], binHashes[x], weight=1)
        
        # Assert
        for connectedBins in nx.connected_components(binGraph):
            connectedBins = list(connectedBins)
            self.assertEqual(len(connectedBins), 3, "Should have 3 connected components")

class TestFragmentMerger(unittest.TestCase):
    def test_fragment_merger(self):
        '''
        This test should result in the fragmented gene bins being merged together
        '''
        # Arrange
        gff3Files = [os.path.join(dataDir, "annotation_fragments.gff3")]
        binCollectionList = _generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0]
        gmapFiles = [os.path.join(dataDir, "gmap_fragments.gff3")]
        
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        
        binCollectionList = _generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0] # offset the change to this test
        
        # Act
        origNumBins = len(binCollection)
        origNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        for i in range(len(binCollectionList)):
            binCollection = binCollectionList[i].fix_fragments(multiOverlaps[i])
            binCollectionList[i] = binCollection
        
        newNumBins = len(binCollection)
        newNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        
        # Assert
        self.assertEqual(origNumBins, 2, "Should contain 2 bins")
        self.assertEqual(newNumBins, 1, "Should contain 1 bin")
        
        self.assertEqual(origNumIDs, [1, 1], "Should contain 2 bins with 1 ID each")
        self.assertEqual(newNumIDs, [4], "Should contain 1 bin with 4 IDs")

class TestBinLinking(unittest.TestCase):
    def test_bin_selfmerger(self):
        '''
        This test should result in the bins (from separate genomes) being merged together.
        '''
        # Arrange
        gmapFiles = [
            os.path.join(dataDir, "gmap_cmj.gff3"),
            os.path.join(dataDir, "gmap_cmj_cross.gff3"),
            os.path.join(dataDir, "gmap_fh.gff3"),
            os.path.join(dataDir, "gmap_fh_cross.gff3")
        ]
        binCollectionList = [ BinCollection() for _ in range(len(gmapFiles)) ]
        # collectionList = binCollectionList
        binCollectionList, multiOverlaps = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1,
                                                                     gmapIdentity=0.95)
        convergenceIters = 5
        
        # Act
        origNumBins = [ len(bc) for bc in binCollectionList ]
        origNumIDs = [ len(bin.data.ids) for bc in binCollectionList for bin in bc ]
        origIDs = [ bin.data.ids for bc in binCollectionList for bin in bc ]
        origExons = [ bin.data.exons for bc in binCollectionList for bin in bc ]
        origContigs = [ bin.data.contig for bc in binCollectionList for bin in bc ]
        
        binCollection = binCollectionList[0]
        for i in range(1, len(binCollectionList)):
            binCollection.merge(binCollectionList[i])
        
        mergedNumBins = len(binCollection)
        mergedNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        mergedIDs = [ bin.data.ids for bin in binCollection ]
        mergedExons = [ bin.data.exons for bin in binCollection ]
        mergedContigs = [ bin.data.contig for bin in binCollection ]
        
        binCollection = iterative_bin_self_linking(binCollection, convergenceIters)
        
        newNumBins = len(binCollection)
        newNumIDs = [ len(bin.data.ids) for bin in binCollection ]
        newIDs = [ bin.data.ids for bin in binCollection ]
        newExons = [ bin.data.exons for bin in binCollection ]
        newContigs = [ bin.data.contig for bin in binCollection ]
        
        # Assert
        self.assertEqual(origNumBins, [1, 1, 1, 1], "Should contain 4 binCollections with 1 bin each")
        self.assertEqual(newNumBins, 2, "Should contain 2 bins")
        self.assertEqual(newNumIDs, [1, 1], "Should contain 2 bins with 1 ID each")

class TestMultiProcessing(unittest.TestCase):
    def test_multiprocessing_equality(self):
        '''
        This test should give the same result regardless of number of threads
        '''
        # Arrange
        convergenceIters = 5
        gmapFiles = [
            os.path.join(dataDir, "gmap_cmj.gff3"),
            os.path.join(dataDir, "gmap_cmj_cross.gff3"),
            os.path.join(dataDir, "gmap_fh.gff3"),
            os.path.join(dataDir, "gmap_fh_cross.gff3")
        ]
        binCollectionList = [ BinCollection() for _ in range(len(gmapFiles)) ]
        
        # Act 1
        binCollectionList1 = [ BinCollection() for _ in range(len(gmapFiles)) ]
        binCollectionList1, multiOverlaps1 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1,
                                                                     gmapIdentity=0.95)
        origIDs1 = [ bin.data.ids for bc in binCollectionList1 for bin in bc ]
        origExons1 = [ bin.data.exons for bc in binCollectionList1 for bin in bc ]
        origContigs1 = [ bin.data.contig for bc in binCollectionList1 for bin in bc ]
        
        # Act 2
        binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
        binCollectionList2, multiOverlaps2 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=4,
                                                                     gmapIdentity=0.95)
        origIDs2 = [ bin.data.ids for bc in binCollectionList2 for bin in bc ]
        origExons2 = [ bin.data.exons for bc in binCollectionList2 for bin in bc ]
        origContigs2 = [ bin.data.contig for bc in binCollectionList2 for bin in bc ]
        
        # Assert
        self.assertEqual(origIDs1, origIDs2, "Thread number should have no effect")
        self.assertEqual(origExons1, origExons2, "Thread number should have no effect")
        self.assertEqual(origContigs1, origContigs2, "Thread number should have no effect")

    def test_multiprocessing_speed(self):
            '''
            This test should prove that multiprocessing has a speed benefit
            '''
            # Arrange
            convergenceIters = 5
            gmapFiles = [
                os.path.join(dataDir, "gmap_cmj.gff3"),
                os.path.join(dataDir, "gmap_cmj_cross.gff3"),
                os.path.join(dataDir, "gmap_fh.gff3"),
                os.path.join(dataDir, "gmap_fh_cross.gff3")
            ]
            binCollectionList = [ BinCollection() for _ in range(len(gmapFiles)) ]
            numTests = 100
            
            # Act 1
            act1Time = 0
            for i in range(numTests):
                timeStart = time.time()
                binCollectionList1 = [ BinCollection() for _ in range(len(gmapFiles)) ]
                binCollectionList1, multiOverlaps1 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                            threads=1,
                                                                            gmapIdentity=0.95)
                act1Time += time.time() - timeStart
            
            # Act 2
            act2Time = 0
            for i in range(numTests):
                timeStart = time.time()
                binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
                binCollectionList2, multiOverlaps2 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                            threads=4,
                                                                            gmapIdentity=0.95)
                act2Time += time.time() - timeStart
            
            # Assert
            self.assertGreater(act1Time, act2Time, "time 2 should be less than time 1")

if __name__ == '__main__':
    unittest.main()
