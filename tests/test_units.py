#!/usr/bin/env python3

import os, sys, unittest
import networkx as nx

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.bins import BinCollection, Bin, BinBundle
from modules.gff3 import GmapGFF3
from modules.bin_handling import GmapBinProcess, CollectionSeedProcess, \
    find_overlapping_bins, add_bin_to_collection
from modules.parsing import BINge_Results

# Specify data locations
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dataDir = os.path.join(baseDir, "tests", "data")

###

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

def generate_bin_collections(gff3Files, threads=1, isMicrobial=False):
    '''
    A re-implementation of the generate_bin_collections function that uses the
    CollectionSeedProcess class to create the BinCollections in parallel.
    
    Parameters:
        gff3Files -- a list of GFF3 files to be processed.
        threads -- an integer indicating how many threads to run.
        isMicrobial -- a boolean indicating whether the genomes are microbial or not which,
                       in turn, determines whether we will parse mRNA features (False) or
                       gene features (True).
    Returns:
        collectionList -- a list containing BinCollections in numerical order of the genome
                          files.
    '''
    # Start up threads
    collectionList = []
    for i in range(0, len(gff3Files), threads): # only process n (threads) at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(gff3Files): # parent loop may excess if n > the number of GMAP files
                gff3File = gff3Files[i+x]
                
                seedWorkerThread = CollectionSeedProcess(gff3File, isMicrobial) # returns empty BinCollection if no GFF3
                seedWorkerThread.start()
                processing.append(seedWorkerThread)
        
        # Wait on processes to end
        for seedWorkerThread in processing:
            binCollection = seedWorkerThread.get_result()
            seedWorkerThread.join()
            seedWorkerThread.check_errors()
            
            collectionList.append(binCollection)
    
    return collectionList

def populate_bin_collections(collectionList, gmapFiles, threads=1):
    '''
    A re-implement with the same logic just minus some of the setup.
    
    Parameters:
        collectionList -- a list of BinCollections to be populated.
        gmapFiles -- a list of GMAP files to be processed.
        threads -- an integer indicating how many threads to run.
    Returns:
        collectionList -- a list of BinCollections, each BinCollection containing
                          Bins populated with GMAP alignments.
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
    collectionList = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection = threadData[i+x]
                indexFile = os.path.join(dataDir, "length_index.pkl")
                
                populateWorkerThread = GmapBinProcess(thisGmapFiles, thisBinCollection,
                                                     indexFile)
                
                processing.append(populateWorkerThread)
                populateWorkerThread.start()
        
        # Gather results
        for populateWorkerThread in processing:
            binCollection = populateWorkerThread.get_result()
            populateWorkerThread.join()
            populateWorkerThread.check_errors()
            collectionList.append(binCollection)
    
    return collectionList

###

# Define unit tests
class TestBinCollection(unittest.TestCase):
    def test_binCollection_init(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "gmap_normal.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0]
        
        # Act
        ## N/A
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should be 2")
    
    def test_binCollection_find(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        binCollection = binCollectionList[0]
        
        # Act
        binOverlap = binCollection.find("genome1", 1, 50)
        binOverlap2 = binCollection.find("genome1", 1, 100)
        
        # Assert
        self.assertEqual(len(binOverlap), 1, "Should be 1")
        
        self.assertEqual(binOverlap[0].start, 1, "Should be 1")
        self.assertEqual(binOverlap[0].end, 50, "Should be 50")
        self.assertEqual(binOverlap[0].contig, "genome1", "Should be genome1")
        
        self.assertEqual(len(binOverlap[0].ids), 1, "Should contain 1")
        self.assertIn("genome1.1.mrna", binOverlap[0].ids, "Should contain genome1.1.mrna")
        
        self.assertEqual(len(binOverlap2), 2, "Should be 2")
    
    def test_binCollection_delete(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList1 = generate_bin_collections(gff3Files)
        binCollectionList2 = generate_bin_collections(gff3Files)
        
        binCollection1 = binCollectionList1[0]
        binCollection2 = binCollectionList2[0]
        
        binOverlap1 = binCollection1.find("genome1", 1, 50)
        binOverlap2 = binCollection2.find("genome1", 51, 100)
        
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
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("genome1", 1, 100)
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
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("genome1", 1, 100)
        bin1 = binOverlap[0]
        bin2 = binOverlap[1]
        
        # Act
        bin1.merge(bin1)
        bin2.merge(bin2)
        
        # Assert
        self.assertEqual(bin2.start, 51, "Should be 51")
        self.assertEqual(bin2.end, 100, "Should be 100")
        self.assertEqual(bin2.contig, "genome1", "Should be genome1")
        
        self.assertEqual(len(bin1.ids), 1, "Should contain 1 id")
        self.assertEqual(len(bin2.ids), 1, "Should contain 1 id")
        
        self.assertIn("genome1.2.mrna", bin2.ids, "Should contain genome1.2.mrna")
    
    def test_bin_merge(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        
        binCollection = binCollectionList[0]
        binOverlap = binCollection.find("genome1", 1, 100)
        bin1 = binOverlap[0]
        bin2 = binOverlap[1]
        
        # Act
        bin1.merge(bin2)
        
        # Assert
        self.assertEqual(bin1.start, 1, "Should be 1")
        self.assertEqual(bin1.end, 100, "Should be 100")
        self.assertEqual(bin1.contig, "genome1", "Should be genome1")
        
        self.assertEqual(len(bin1.ids), 2, "Should contain 2 ids")
        
        self.assertIn("genome1.1.mrna", bin1.ids, "Should contain genome1.1.mrna")   
        self.assertIn("genome1.2.mrna", bin1.ids, "Should contain genome1.2.mrna")

class TestGff3Iterate(unittest.TestCase):
    def test_iterate_through_empty_file(self):
        # Arrange
        gmapFile = os.path.join(dataDir, "empty_file")
        gmapGff3Obj = GmapGFF3(gmapFile)
        
        # Act and assert
        try:
            for feature in gmapGff3Obj:
                pass
        except:
            self.assertTrue(True, "Should error out here")
    
    def test_iterate_through_gff3(self):
        # Arrange
        gmapFile = os.path.join(dataDir, "gmap_normal.gff3")
        gmapGff3Obj = GmapGFF3(gmapFile)
        
        # Act
        numFeatures = 0
        mrnaIDs = []
        for feature in gmapGff3Obj:
            numFeatures += 1
            mrnaIDs.append(feature["Name"])
        
        # Assert
        self.assertEqual(numFeatures, 2, "Should contain 2 bins")
        self.assertEqual(mrnaIDs, ['genome1.1', 'genome1.2'],
                         "Should be ['genome1.1', 'genome1.2']")

class TestNovelPopulate(unittest.TestCase):
    def test_populate_reject_novel(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        binCollectionList = populate_bin_collections(binCollectionList, gmapFiles,
                                                     threads=1)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 0, "Should contain 0 bins")
    
    def test_populate_annotation_only(self):
        # Arrange
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        binCollectionList = populate_bin_collections(binCollectionList, gmapFiles,
                                                     threads=1)
        binCollection = binCollectionList[0]
        
        binIDs = [ x for bin in binCollection.bins["genome1"] for x in bin.data.ids ]
        
        # Assert
        self.assertEqual(len(binCollection), 9, "Should contain 9 bins")
        self.assertTrue(not any([x.startswith("contig1") for x in binIDs]),
                        "Should not contain any contig1 IDs")
    
    def test_populate_novel_only(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollectionList = populate_bin_collections(binCollectionList, gmapFiles,
                                                     threads=1)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should contain 2 bins")
    
    def test_populate_novel_merge(self):
        "I think this test isn't really testing what it's meant to anymore"
        # Arrange
        gff3Files = [os.path.join(dataDir, "genome1.gff3")]
        binCollectionList = generate_bin_collections(gff3Files)
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        binCollectionList = populate_bin_collections(binCollectionList, gmapFiles,
                                                     threads=1)
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 9, "Should contain 6 bins")
        for bin in binCollection:
            if bin.data.start == 1 or bin.data.start == 51:
                self.assertEqual(len(bin.data.ids), 2, f"Should have 2 IDs in bin")

class TestNetworkX(unittest.TestCase):
    def test_connected_component_numbers(self):
        # Arrange
        binCollection = BinCollection()
        bin1 = Bin("genome1", 1, 100)
        bin1.ids = {"genome1.1"}
        bin1.exons = {"genome1.1": [[1, 100]]}
        
        bin2 = Bin("genome1", 99, 200)
        bin2.ids = {"genome1.2"}
        bin2.exons = {"genome1.2": [[99, 100200]]}
        
        bin3 = Bin("genome1", 199, 300)
        bin3.ids = {"genome1.3"}
        bin3.exons = {"genome1.3": [[199, 300]]}
        
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
        binCollectionList1 = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1)
        origIDs1 = [ bin.data.ids for bc in binCollectionList1 for bin in bc ]
        origStart1 = [ bin.data.start for bc in binCollectionList1 for bin in bc ]
        origContigs1 = [ bin.data.contig for bc in binCollectionList1 for bin in bc ]
        
        # Act 2
        binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
        binCollectionList2 = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=4)
        origIDs2 = [ bin.data.ids for bc in binCollectionList2 for bin in bc ]
        origStart2 = [ bin.data.start for bc in binCollectionList2 for bin in bc ]
        origContigs2 = [ bin.data.contig for bc in binCollectionList2 for bin in bc ]
        
        # Assert
        self.assertEqual(origIDs1, origIDs2, "Thread number should have no effect")
        self.assertEqual(origStart1, origStart2, "Thread number should have no effect")
        self.assertEqual(origContigs1, origContigs2, "Thread number should have no effect")
    
    ## Turn off this test in case user has a single core machine; it will probably fail
    # def test_multiprocessing_speed(self):
    #         '''
    #         This test should prove that multiprocessing has a speed benefit
    #         '''
    #         # Arrange
    #         convergenceIters = 5
    #         gmapFiles = [
    #             os.path.join(dataDir, "gmap_cmj.gff3"),
    #             os.path.join(dataDir, "gmap_cmj_cross.gff3"),
    #             os.path.join(dataDir, "gmap_fh.gff3"),
    #             os.path.join(dataDir, "gmap_fh_cross.gff3")
    #         ]
    #         binCollectionList = [ BinCollection() for _ in range(len(gmapFiles)) ]
    #         numTests = 100
            
    #         # Act 1
    #         act1Time = 0
    #         for i in range(numTests):
    #             timeStart = time.time()
    #             binCollectionList1 = [ BinCollection() for _ in range(len(gmapFiles)) ]
    #             binCollectionList1 = populate_bin_collections(binCollectionList, gmapFiles,
    #                                                                         threads=1)
    #             act1Time += time.time() - timeStart
            
    #         # Act 2
    #         act2Time = 0
    #         for i in range(numTests):
    #             timeStart = time.time()
    #             binCollectionList2 = [ BinCollection() for _ in range(len(gmapFiles)) ]
    #             binCollectionList2 = populate_bin_collections(binCollectionList, gmapFiles,
    #                                                                         threads=4)
    #             act2Time += time.time() - timeStart
            
    #         # Assert
    #         self.assertGreater(act1Time, act2Time, "time 2 should be less than time 1")

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
        indexFile = os.path.join(dataDir, "length_index.pkl")
        binCollection = BinCollection()
        
        # Act
        processor = GmapBinProcess(gff3Files, binCollection, indexFile)
        #gmapFiles, binCollection, indexFileName = gff3Files, binCollection, indexFile
        processor.start()
        resultCollection = processor.get_result()
        processor.join()
        processor.check_errors()
        
        # Assert
        self.assertEqual(len(resultCollection), 2, "Should have two bins")

class TestBINge_Results(unittest.TestCase):
    def test_parse(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        trueBinnedIDs = [f"g{i}" for i in range(1, 8)] # g1 through g7
        trueUnbinnedIDs = [f"g{i}" for i in range(8, 11)] # g8 through g10
        
        # Act
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        binnedIDs = [ v for values in bingeResults.binned.values() for v in values]
        unbinnedIDs = [ v for values in bingeResults.unbinned.values() for v in values]
        
        # Assert
        self.assertEqual(len(bingeResults.binned), 5, "Should have five binned clusters")
        self.assertEqual(len(bingeResults.unbinned), 2, "Should have two binned clusters")
        self.assertEqual(binnedIDs, trueBinnedIDs, f"Binned IDs should be {trueBinnedIDs}, not {binnedIDs}")
        self.assertEqual(unbinnedIDs, trueUnbinnedIDs, f"Unbinned IDs should be {trueUnbinnedIDs}, not {unbinnedIDs}")
    
    def test_set_1(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        
        # Act
        bingeResults2 = BINge_Results()
        bingeResults2.binned = bingeResults.binned
        bingeResults2.unbinned = bingeResults.unbinned
        
        # Assert
        self.assertEqual(bingeResults.binned, bingeResults2.binned, "Binned .parse should be equal to set data")
        self.assertEqual(bingeResults.unbinned, bingeResults2.unbinned, "Unbinned .parse should be equal to set data")
    
    def test_set_2(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        
        binnedClusters = {0: ['g1'], 1: ['g2'], 2: ['g3'], 3: ['g4', 'g5'], 4: ['g6', 'g7']}
        unbinnedClusters = {0: ['g8', 'g9'], 1: ['g10']}
        
        # Act
        bingeResults2 = BINge_Results()
        bingeResults2.binned = binnedClusters
        bingeResults2.unbinned = bingeResults2.update_unbinned_ids(unbinnedClusters)
        
        # Assert
        self.assertEqual(bingeResults.binned, bingeResults2.binned, "Binned .parse should be equal to set data")
        self.assertEqual(bingeResults.unbinned, bingeResults2.unbinned, "Unbinned .parse should be equal to set data")
    
    def test_filter(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        
        binnedClusterKeys = [1,2,3,4,5]
        unbinnedClusterKeys = [6]
        toRemove = [0, 6] # 0: ['g1'] and 6: ['g10'] should be removed
        
        # Act
        bingeResults.filter(toRemove)
        
        # Assert
        self.assertEqual(len(bingeResults.binned), 4, "Should have four binned clusters")
        self.assertEqual(len(bingeResults.unbinned), 1, "Should have one binned cluster")
    
    def test_iterate(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        
        trueClusterIDs = list(range(0, 7)) # [0,1,2,3,4,5,6,]
        trueSeqIDs = [f"g{i}" for i in range(1, 11)] # g1 through g10
        
        # Act
        foundClusterIDs = []
        foundSeqIDs = []
        for clusterID, seqIDs in bingeResults:
            foundClusterIDs.append(int(clusterID))
            foundSeqIDs.extend(seqIDs)
        
        # Assert
        self.assertEqual(foundClusterIDs, trueClusterIDs, f"Cluster IDs should be {trueClusterIDs}, not {foundClusterIDs}")
        self.assertEqual(foundSeqIDs, trueSeqIDs, f"Sequence IDs should be {trueSeqIDs}, not {foundSeqIDs}")
    
    def test_write(self):
        # Arrange
        resultFile = os.path.join(dataDir, "clustering_result_1.tsv")
        bingeResults = BINge_Results()
        bingeResults.parse(resultFile)
        
        tmpFile = os.path.join(dataDir, "TestBINge_Results.test_write.tsv")
        
        # Act
        bingeResults.write(tmpFile)
        bingeResults2 = BINge_Results()
        bingeResults2.parse(tmpFile)
        
        # Assert
        self.assertEqual(bingeResults.binned, bingeResults2.binned, "Binned .write should be equal to .parse data")
        self.assertEqual(bingeResults.unbinned, bingeResults2.unbinned, "Unbinned .write should be equal to .parse data")

if __name__ == '__main__':
    unittest.main()
