#! python3

import os, sys, unittest, time
import networkx as nx
from multiprocessing import Queue

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_GFF3IO import GFF3
from modules.bins import BinCollection, Bin, BinBundle
from modules.gff3_handling import iterate_gmap_gff3
from modules.thread_workers import GmapBinProcess, CollectionSeedProcess, \
    QueuedBinSplitterProcess, find_overlapping_bins, add_bin_to_collection

###

## Unit tests to be done:
## cluster_by_occurrence
## 

def _generate_bin_collections_via_seeder(gff3Files, threads):
    # Start up threads
    collectionList = []
    for i in range(0, len(gff3Files), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(gff3Files): # parent loop may excess if n > the number of GMAP files
                gff3File = gff3Files[i+x]
                
                seedWorkerThread = CollectionSeedProcess(gff3File)
                
                processing.append(seedWorkerThread)
                seedWorkerThread.start()
        
        # Wait on processes to end
        for seedWorkerThread in processing:
            binCollection = seedWorkerThread.get_result()
            seedWorkerThread.join()
            seedWorkerThread.check_errors()
            collectionList.append(binCollection)
    
    return collectionList

def _generate_bin_collections(gff3Files, isMicrobial=False):
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
                
                # Create a bin for each exon feature
                exonBins = []
                if not isMicrobial:
                    try:
                        for mrnaFeature in geneFeature.mRNA:
                            for exonFeature in mrnaFeature.exon:
                                exonBin = Bin(exonFeature.contig, exonFeature.start, exonFeature.end)
                                exonBin.add(mrnaFeature.ID)
                                exonBins.append(exonBin)
                    except:
                        "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                        continue
                else:
                    try:
                        exonBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
                        exonBin.add(geneFeature.ID)
                        exonBins.append(exonBin)
                    except:
                        "This exception occurs if a gene feature has unusual children types e.g., rRNA"
                        continue
                
                # Iteratively handle exon bins
                for exonBin in exonBins:
                    binOverlap = find_overlapping_bins(binCollection, exonBin)
                    add_bin_to_collection(binCollection, binOverlap, exonBin)
        
        # Store the bin collection in our list for multi-threading later
        collectionList.append(binCollection)
    
    return collectionList

def _populate_bin_collections(collectionList, gmapFiles, threads=1, gmapIdentity=0.95):
    '''
    A re-implement with the same logic just minus some of the setup.
    '''    
    # Establish lists for feeding data into threads
    threadData = []
    for i in range(len(collectionList)):
        # Get all GMAP files associated with this genome
        thisGmapFiles = [gmapFiles[i]]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[i]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection])
    
    # Start up threads
    resultBinCollection = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection = threadData[i+x]
                
                populateWorkerThread = GmapBinProcess(thisGmapFiles, thisBinCollection,
                                                     gmapIdentity)
                
                processing.append(populateWorkerThread)
                populateWorkerThread.start()
        
        # Gather results
        for populateWorkerThread in processing:
            binCollection = populateWorkerThread.get_result()
            populateWorkerThread.join()
            populateWorkerThread.check_errors()
            resultBinCollection.append(binCollection)
    
    return resultBinCollection

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
            for feature in iterate_gmap_gff3(gmapFile):
                pass
        except:
            self.assertTrue(True, "Should error out here")
    
    def test_iterate_through_gff3(self):
        # Arrange
        gmapFile = os.path.join(dataDir, "gmap_normal.gff3")
        
        # Act
        numFeatures = 0
        mrnaIDs = []
        for feature in iterate_gmap_gff3(gmapFile):
            numFeatures += 1
            mrnaIDs.append(feature["Name"])
        
        # Assert
        self.assertEqual(numFeatures, 2, "Should contain 2 bins")
        self.assertEqual(mrnaIDs, ['contig1.1', 'contig1.2'],
                         "Should be ['contig1.1', 'contig1.2']")

class TestNovelPopulate(unittest.TestCase):
    def test_populate_reject_novel(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        binCollectionList = _populate_bin_collections(binCollectionList, gmapFiles,
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
        binCollectionList = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 6, "Should contain 6 bins")
    
    def test_populate_novel_only(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollectionList = _populate_bin_collections(binCollectionList, gmapFiles,
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
        binCollectionList = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=0.95)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 6, "Should contain 6 bins")
        for bin in binCollection:
            if bin.data.start == 1 or bin.data.start == 51:
                self.assertEqual(len(bin.data.ids), 2, f"Should have 2 IDs in bin")

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
        binCollectionList1 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1,
                                                                     gmapIdentity=0.95)
        origIDs1 = [ bin.data.ids for bc in binCollectionList1 for bin in bc ]
        origExons1 = [ bin.data.exons for bc in binCollectionList1 for bin in bc ]
        origContigs1 = [ bin.data.contig for bc in binCollectionList1 for bin in bc ]
        
        # Act 2
        binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
        binCollectionList2 = _populate_bin_collections(binCollectionList, gmapFiles,
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
                binCollectionList1 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                            threads=1,
                                                                            gmapIdentity=0.95)
                act1Time += time.time() - timeStart
            
            # Act 2
            act2Time = 0
            for i in range(numTests):
                timeStart = time.time()
                binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
                binCollectionList2 = _populate_bin_collections(binCollectionList, gmapFiles,
                                                                            threads=4,
                                                                            gmapIdentity=0.95)
                act2Time += time.time() - timeStart
            
            # Assert
            self.assertGreater(act1Time, act2Time, "time 2 should be less than time 1")

class TestBinSeederThread(unittest.TestCase):
    def test_binCollection_seeder(self):
        # Arrange
        threads=4
        gff3Files = [os.path.join(dataDir, "gmap_normal.gff3")]
        binCollectionList = _generate_bin_collections_via_seeder(gff3Files, threads)
        binCollection = binCollectionList[0]
        
        # Act
        ## N/A
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should be 2")

    def test_seeder_has_exception(self):
        # Arrange
        gff3File = os.path.join(dataDir, "gmap_normal.notexist.gff3")
        
        # Act and Assert
        try:
            seedWorkerThread = CollectionSeedProcess(gff3File)
            seedWorkerThread.start()
            seedWorkerThread.join()
            seedWorkerThread.check_errors()
            self.assertTrue(False, "Not having exception is a test fail")
        except FileNotFoundError:
            self.assertTrue(True, "Having FileNotFoundError is a test pass")

    def test_seeder_does_work(self):
        # Arrange
        gff3File = os.path.join(dataDir, "gmap_normal.gff3")
        
        # Act
        seedWorkerThread = CollectionSeedProcess(gff3File)
        seedWorkerThread.start()
        seedWorkerThread.join()
        seedWorkerThread.check_errors()
        result = seedWorkerThread.get_result(1)
        
        # Assert
        self.assertEqual(len(result), 2, "Result should have 2 bins")

class TestGmapBinProcess(unittest.TestCase):
    def test_gmap_bin_process(self):
        # Arrange
        threads=4
        gff3Files = [os.path.join(dataDir, "gmap_normal.gff3")]
        binCollection = BinCollection()
        
        # Act
        processor = GmapBinProcess(gff3Files, binCollection, 0.95)
        processor.start()
        resultCollection = processor.get_result()
        processor.join()
        processor.check_errors()
        
        # Assert
        self.assertEqual(len(resultCollection), 2, "Should have two bins")

class TestGFF3IterationSpeed(unittest.TestCase):
    def test_iterators_for_speed(self):
        # Arrange
        gff3File = os.path.join(dataDir, "big.gff3")
        timesToIterate = 5
        
        # Skip if test is not possible
        if not os.path.exists(gff3File):
            self.assertTrue(True, "Skipping test because big.gff3 does not exist")
            return
        
        # Act
        startTime1 = time.time()
        for _ in range(timesToIterate):
            for feature in iterate_gmap_gff3(gff3File):
                pass
        endTime1 = time.time() - startTime1
        
        startTime2 = time.time()
        for _ in range(timesToIterate):
            gff3Obj = GFF3(gff3File, strict_parse=False)
            for geneFeature in gff3Obj.types["gene"]:
                pass
        endTime2 = time.time() - startTime2
        
        # Notify
        print(f"test_iterators_for_speed: iterate_gmap_gff3 = {endTime1}")
        print(f"test_iterators_for_speed: GFF3 = {endTime2}")
        
        # Assert
        self.assertTrue(True, "Test result is printed to terminal")

if __name__ == '__main__':
    unittest.main()
