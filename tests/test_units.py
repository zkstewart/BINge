#! python3

import os, sys, unittest
import networkx as nx

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_GFF3IO
from modules.fasta_handling import FastaCollection
from modules.bins import BinCollection, Bin
from modules.gff3_handling import iterate_through_gff3
from BINge import generate_bin_collections, populate_bin_collections, \
    multithread_bin_splitter, iterative_bin_self_linking
from test_scenarios import get_binCollection_in_range, binge_runner

# Specify data locations
dataDir = os.path.join(os.getcwd(), "data")

# Define unit tests
class TestBinCollection(unittest.TestCase):
    def test_binCollection_init(self):
        # Arrange
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "gmap_normal.gff3")], 1)
        binCollection = binCollectionList[0]
        
        # Act
        ## N/A
        
        # Assert
        self.assertEqual(len(binCollection), 2, "Should be 2")
    
    def test_binCollection_find(self):
        # Arrange
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
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
        self.assertIn("contig1.1", binOverlap[0].ids, "Should contain contig1.1")
        
        self.assertEqual(len(binOverlap[0].exons), 1, "Should contain 1")
        self.assertIn("contig1.1", binOverlap[0].exons, "Should contain contig1.1")
        self.assertEqual(binOverlap[0].exons["contig1.1"], [[1, 50]], "Should be [[1, 50]]")
        
        self.assertEqual(len(binOverlap2), 2, "Should be 2")
    
    def test_binCollection_delete(self):
        # Arrange
        binCollectionList1 = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        binCollectionList2 = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        
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
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        
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
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        
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
        
        self.assertIn("contig1.2", bin2.ids, "Should contain contig1.2")
    
    def test_bin_merge(self):
        # Arrange
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        
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
        
        self.assertIn("contig1.1", bin1.ids, "Should contain contig1.1")   
        self.assertIn("contig1.2", bin1.ids, "Should contain contig1.2")

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
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=[0.95])
        
        # Assert
        self.assertEqual(len(novelBinCollection), 0, "Should contain 0 bins")
    
    def test_populate_annotation_only(self):
        # Arrange
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        gmapFiles = [os.path.join(dataDir, "gmap_bad.gff3")]
        
        # Act
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=[0.95])
        binCollection = binCollectionList[0]
        
        # Assert
        self.assertEqual(len(binCollection), 6, "Should contain 6 bins")
    
    def test_populate_novel_only(self):
        # Arrange
        binCollectionList = [BinCollection()]
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=[0.95])
        
        # Assert
        self.assertEqual(len(novelBinCollection), 2, "Should contain 2 bins")
    
    def test_populate_novel_merge(self):
        # Arrange
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation.gff3")], 1)
        gmapFiles = [os.path.join(dataDir, "gmap_normal.gff3")]
        
        # Act
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=[0.95])
        
        # Assert
        self.assertEqual(len(novelBinCollection), 0, "Should contain 0 bins")
        for bin in binCollectionList[0]:
            if bin.data.start == 1 or bin.data.start == 51:
                self.assertEqual(len(bin.data.ids), 1, f"Should have 1 ID in bin")

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
            binCollection, novelBinCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
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
        binCollection, novelBinCollection = binge_runner(binCollectionList, gmapFiles, threads=1, convergenceIters=5)
        
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
        binCollectionList = generate_bin_collections([os.path.join(dataDir, "annotation_fragments.gff3")], 1)
        binCollection = binCollectionList[0]
        gmapFiles = [os.path.join(dataDir, "gmap_fragments.gff3")]
        
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1, gmapIdentity=[0.95])
        
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
        
        novelBinCollection, multiOverlaps = populate_bin_collections(binCollectionList, gmapFiles,
                                                                     threads=1,
                                                                     gmapIdentity=[0.95]*len(gmapFiles))
        convergenceIters = 5
        
        # Act
        origNumBins = len(novelBinCollection)
        origNumIDs = [ len(bin.data.ids) for bin in novelBinCollection ]
        #origIDs = [ bin.data.ids for bin in novelBinCollection ]
        #origExons = [ bin.data.exons for bin in novelBinCollection ]
        #origContigs = [ bin.data.contig for bin in novelBinCollection ]
        
        novelBinCollection = iterative_bin_self_linking(novelBinCollection, convergenceIters)
        
        newNumBins = len(novelBinCollection)
        newNumIDs = [ len(bin.data.ids) for bin in novelBinCollection ]
        #newIDs = [ bin.data.ids for bin in novelBinCollection ]
        #newExons = [ bin.data.exons for bin in novelBinCollection ]
        #newContigs = [ bin.data.contig for bin in novelBinCollection ]
        
        # Assert
        self.assertEqual(origNumBins, 2, "Should contain 2 bins")
        self.assertEqual(newNumBins, 1, "Should contain 1 bin")
        self.assertEqual(origNumIDs, [2, 2], "Should contain 2 bins with 2 IDs each")
        self.assertEqual(newNumIDs, [2], "Should contain 1 bin with 2 IDs")

if __name__ == '__main__':
    unittest.main()
