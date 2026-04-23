#!/usr/bin/env python3

import os, sys, shutil, subprocess, unittest
import networkx as nx
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.bins import BinCollection, Bin, BinBundle
from modules.gff3 import GmapGFF3, GFF3Graph
from modules.gff3tofasta import gff3_to_fasta
from modules.fasta_handling import txome_to_orfs
from modules.bin_handling import GmapBinProcess, CollectionSeedProcess, \
    find_overlapping_bins, add_bin_to_collection
from modules.gmap_handling import GMAP
from modules.parsing import BINge_Results, load_sequence_length_index
from modules.setup import TargetGenome, AnnotatedGenome, Transcriptome, \
    inputs_to_json, json_to_inputs
from modules.locations import Locations
from modules.validation import check_for_duplicates, check_for_seqid_consistency

# Specify data locations
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dataDir = os.path.join(baseDir, "tests", "data")
workDir = os.path.join(dataDir, "tmp")

###

# Define main program behaviour reimplementations
"these reimplementations should be functionally equivalent, just easier to work with in a test suite"
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

# Define utility functions to enable testing
def run_subprocess(command):
    '''
    Parameters:
        command -- a list of strings representing the command to run
    '''
    process = subprocess.Popen(" ".join(command), shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()
    return process.returncode, stdout.decode("utf-8").rstrip("\r\n "), stderr.decode("utf-8").rstrip("\r\n ")

def cleanup_working_directory():
    # Clean any existing workdir and setup a fresh dir
    if os.path.exists(workDir):
        shutil.rmtree(workDir)
    if not os.path.exists(workDir):
        os.makedirs(workDir)

def make_working_directory(locations):
    os.makedirs(locations.sequencesDir, exist_ok=True)
    os.makedirs(locations.gff3Dir, exist_ok=True)
    os.makedirs(locations.txDir, exist_ok=True)
    os.makedirs(locations.genomesDir, exist_ok=True)
    os.makedirs(locations.mappingDir, exist_ok=True)
    os.makedirs(locations.analysisDir, exist_ok=True)

def obtain_input_args1():
    i_gff1 = os.path.join(dataDir, "genome1.gff3")
    i_fasta1 = os.path.join(dataDir, "genome1.fasta")
    
    i_gff2 = os.path.join(dataDir, "genome2.gff3")
    i_fasta2 = os.path.join(dataDir, "genome2.fasta")
    
    ig_gff1 = os.path.join(dataDir, "genome1.gff3")
    ig_fasta1 = os.path.join(dataDir, "genome1.fasta")
    
    ig_gff2 = os.path.join(dataDir, "genome2.gff3")
    ig_fasta2 = os.path.join(dataDir, "genome2.fasta")
    
    ix_mrna = os.path.join(dataDir, "genome3.ttable1.mrna")
    ix_cds = os.path.join(dataDir, "genome3.ttable1.cds")
    ix_aa = os.path.join(dataDir, "genome3.ttable1.aa")
    
    return i_gff1, i_fasta1, i_gff2, i_fasta2, ig_gff1, ig_fasta1, ig_gff2, ig_fasta2, ix_mrna, ix_cds, ix_aa

def obtain_input_args2():
    'Alternate input arguments which allow (e.g.) --ix args to be extracted without errors'
    gff4 = os.path.join(dataDir, "genome4.gff3")
    fasta4 = os.path.join(dataDir, "genome4.fasta")
    
    gff5 = os.path.join(dataDir, "genome5.gff3")
    fasta5 = os.path.join(dataDir, "genome5.fasta")
    
    mrna6 = os.path.join(dataDir, "genome6.ttable1.mrna")
    cds6 = os.path.join(dataDir, "genome6.ttable1.cds")
    aa6 = os.path.join(dataDir, "genome6.ttable1.aa")
    
    return gff4, fasta4, gff5, fasta5, mrna6, cds6, aa6

def setup_working_directory(targetHasGFF3=True, txomeIsPreExtracted=True, translationTable=1):
    '''
    Replicates the behaviour of the BINge.py's equivalently named function, but
    tailored to this testing suite.
    '''
    cleanup_working_directory()
    
    # Establish locations variable and its downstream folders
    locations = Locations(workDir)
    make_working_directory(locations)
    
    # Establish variables
    i_gff1, i_fasta1, i_gff2, i_fasta2, ig_gff1, ig_fasta1, ig_gff2, ig_fasta2, ix_mrna, ix_cds, ix_aa = \
        obtain_input_args1()
    
    targetGenomes = [
        TargetGenome("1", locations, i_fasta1,
                     gff3=i_gff1 if targetHasGFF3 else None),
        TargetGenome(2, locations, i_fasta2,
                     gff3=i_gff2 if targetHasGFF3 else None)
    ]
    annotatedGenomes = [
        AnnotatedGenome(1, locations, ig_fasta1, ig_gff1),
        AnnotatedGenome("2", locations, ig_fasta2, ig_gff2)
    ]
    transcriptomes = [
        Transcriptome("1", locations, mrnaFasta=ix_mrna,
                                      cdsFasta=ix_cds if txomeIsPreExtracted else None,
                                      protFasta=ix_aa if txomeIsPreExtracted else None)
    ]
    
    # Extract sequences
    for targetGenome in targetGenomes:
        targetGenome.extract_sequences(False, translationTable) # isMicrobial==False
    for annotatedGenome in annotatedGenomes:
        annotatedGenome.extract_sequences(False, translationTable) # isMicrobial==False
    transcriptomes[0].extract_sequences(1) # translationTable==1; need to extract to set .ok flags
    
    return targetGenomes, annotatedGenomes, transcriptomes, locations

###

# Define long unit tests
## Comment these out if undergoing development
class TestGmap(unittest.TestCase):
    def test_gmap(self):
        '''
        Test both TargetGenome index and mapping to prevent longer test times
        '''
        # Arrange #1
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        gmap = shutil.which("gmap")
        gmap_build = shutil.which("gmap_build")
        if gmap == None or gmap_build == None:
            return # do not test if gmap is not in the PATH variable
        gmapDir = os.path.dirname(gmap)
        
        # Act #1
        targetGenome = targetGenomes[0] # just the first one
        targetGenome.gmap_index(gmapDir)
        
        # Assert #1
        self.assertTrue(os.path.isdir(targetGenome.gmap_db), "gmap_build seems to have failed")
        self.assertTrue(os.path.isfile(targetGenome.gmap_db + ".ok"), "gmap_build seems to have failed")
        
        # Arrange #2
        queryFile = os.path.join(dataDir, "genome3.ttable1.mrna")
        outputFileName = os.path.join(locations.mappingDir, "test1.gmap.gff3")
        
        # Act #2
        gmapper = GMAP(queryFile, targetGenome.gmap_db, gmapDir, 1) # threads == 1
        gmapper.gmap(outputFileName)
        
        # Assert #2
        self.assertTrue(os.path.isfile(outputFileName), f"gmap seems to have failed to align '{queryFile}' to '{targetGenome.gmap_db}'")

class TestMain(unittest.TestCase):
    '''
    The tests under this class are intended to run BINge through full system exercises
    rather than as discrete units.
    '''
    def test_with_ipair_noig_noix(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        i_gff1, i_fasta1, _, _, _, _, _, _, _, _, _ = \
            obtain_input_args1()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap")
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{i_gff1},{i_fasta1}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        
        # Cluster
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Cluster assert
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
    
    def test_with_ipair_ig_noix(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        i_gff1, i_fasta1, _, _, ig_gff1, ig_fasta1, ig_gff2, ig_fasta2, _, _, _ = \
            obtain_input_args1()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap")
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{i_gff1},{i_fasta1}",
            "--ig", f"{ig_gff2},{ig_fasta2}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "annotations1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'annotations1' to 'genome1' to have a .ok file")
        
        # Cluster
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Cluster assert
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
    
    def test_with_ipair_ig_ix_extract(self):
        '''This test can't be implemented since we end up trying to use MMseqs2 to cluster
        two unbinned sequences (couldn't map because low complexity) which are shorter
        than 14 bp in length which triggers a MMseqs2 error which cannot be handled
        by BINge'''
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        gff4, fasta4, gff5, fasta5, mrna6, _, _ = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{gff4},{fasta4}",
            "--ig", f"{gff5},{fasta5}",
            "--ix", mrna6,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_ipair_ig_ix_symlink(self):
        '''Same scenario as test_with_ipair_ig_ix_extract()'''
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        gff4, fasta4, gff5, fasta5, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{gff4},{fasta4}",
            "--ig", f"{gff5},{fasta5}",
            "--ix", f"{mrna6},{cds6},{aa6}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_ipair_noig_ix_extract_1(self):
        '''
        This test was originally intended to check normal behaviour, but now with the empty file error
        checking in, it can remain useful as a full system test that it gets triggered during init as
        expected.
        '''
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        i_gff1, i_fasta1, _, _, _, _, _, _, ix_mrna, _, _ = \
            obtain_input_args1()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap")
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{i_gff1},{i_fasta1}",
            "--ix", ix_mrna, # ix_mrna can't be extracted as the ORF is too short and has no ATG start codon
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertTrue("is the file empty or does it not contain any valid open reading frames?" in stderr,
                        f"Should have been unable to extract an ORF from '{ix_mrna}'")
    
    def test_with_ipair_noig_ix_extract_2(self):
        '''
        This test is to perform the original purpose of test_with_ipair_noig_ix_extract_1()
        i.e., testing of --ix when it should extract a gene sequence.
        '''        
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        gff4, fasta4, gff5, fasta5, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{gff4},{fasta4}",
            "--ix", mrna6,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome4' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_ipair_noig_ix_symlink(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        gff4, fasta4, gff5, fasta5, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{gff4},{fasta4}",
            "--ix", f"{mrna6},{cds6},{aa6}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "genome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP GFF3 for 'genome1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_isingle_ig_noix(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        _, i_fasta1, _, _, ig_gff1, ig_fasta1, ig_gff2, ig_fasta2, _, _, _ = \
            obtain_input_args1()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap")
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", i_fasta1,
            "--ig", f"{ig_gff2},{ig_fasta2}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "annotations1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'annotations1' to 'genome1' to have a .ok file")
        
        # Cluster
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Cluster assert
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
    
    def test_with_isingle_ig_ix_extract(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        _, fasta4, gff5, fasta5, mrna6, _, _ = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", fasta4,
            "--ig", f"{gff5},{fasta5}",
            "--ix", mrna6,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "annotations1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'annotations1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_isingle_ig_ix_symlink(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        _, fasta4, gff5, fasta5, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", fasta4,
            "--ig", f"{gff5},{fasta5}",
            "--ix", f"{mrna6},{cds6},{aa6}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "annotations1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'annotations1' to 'genome1' to have a .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_isingle_noig_ix_extract(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        _, fasta4, _, _, mrna6, _, _ = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", fasta4,
            "--ix", mrna6,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_with_isingle_noig_ix_symlink(self):
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        
        _, fasta4, _, _, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", fasta4,
            "--ix", f"{mrna6},{cds6},{aa6}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(fasta1_gmap_db + ".ok"),
                        f"Expected GMAP dir for 'genome1' to have an .ok file")
        self.assertTrue(os.path.isfile(os.path.join(locations.mappingDir, "transcriptome1_to_genome1_gmap.gff3.ok")),
                        f"Expected GMAP mapping for 'transcriptome1' to 'genome1' to have a .ok file")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        clusterFile = os.path.join(locations.analysisDir, locations.runName, locations.clusterFile)
        bingeResults = BINge_Results()
        bingeResults.parse(clusterFile)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(clusterFile + ".ok"),
                        f"Expected '{clusterFile}' to have an .ok file")
        self.assertTrue(len(bingeResults.unbinned) == 0, "All sequences should have been binned")
        self.assertTrue(len(bingeResults.binned) == 1, "All sequences should have been binned into one cluster")
    
    def test_full_1(self):
        'Full system test of all BINge and BINge_post modules'
        # Init arrange
        cleanup_working_directory()
        locations = Locations(workDir)
        
        gff4, fasta4, gff5, fasta5, mrna6, cds6, aa6 = obtain_input_args2()
        fasta1_gmap_db = os.path.join(locations.genomesDir, "genome1.gmap") # note: it's still called genome1
        threads = 1
        
        # Init act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "initialise",
            "-d", workDir,
            "-i", f"{gff4},{fasta4}",
            "--ix", f"{mrna6},{cds6},{aa6}",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Init assert
        self.assertTrue(stderr == "", f"'init' has stderr output: {stderr}")
        
        # Cluster act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge.py"), "cluster",
            "-d", workDir,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Cluster assert
        self.assertTrue(stderr == "", f"'cluster' has stderr output: {stderr}")
        
        # BLAST arrange
        queryFile = os.path.join(dataDir, "example_uniref90.fasta")
        blastFile = os.path.join(locations.blastDir, locations.blastFile)
        
        # BLAST act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "blast",
            "-d", workDir,
            "-t", queryFile,
            "-s", "protein",
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # BLAST assert
        self.assertTrue(stderr == "", f"'blast' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(blastFile + ".ok"),
                        f"Expected '{blastFile}' to have an .ok file")
        
        # Salmon arrange
        readSuffix = "P.fastq"
        salmonResultFile = os.path.join(locations.salmonDir, "genome4_rnaseq")
        quantFile = os.path.join(salmonResultFile, "quant.sf")
        
        # Salmon act
        ## Test reads must have at least 10 sequences or --minAssignedFrags error occurs
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "salmon",
            "-d", workDir,
            "-r", dataDir,
            "-s", readSuffix,
            "--threads", str(threads)
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Salmon assert
        self.assertTrue(stderr == "", f"'salmon' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(salmonResultFile + ".ok"),
                        f"Expected '{salmonResultFile}' to have an .ok file")
        
        with open(quantFile, "r") as fileIn:
            quantContents = fileIn.read()
        numReads = float(quantContents.rstrip().split("\t")[-1])
        
        self.assertTrue(numReads == 10.0, f"quant.sf should have 10.0 reads, not {numReads}")
        
        # Filter arrange
        filterFile = os.path.join(locations.filterDir, locations.runName, locations.filteredClusterFile)
        
        # Filter act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "filter",
            "-d", workDir,
            "--useGFF3", "--useBLAST", "--useSalmon",
            "--readLength", "100"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Filter assert
        self.assertTrue(stderr == "", f"'filter' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(filterFile + ".ok"),
                        f"Expected '{filterFile}' to have an .ok file")
        
        # Representatives arrange
        repMRNAFile = os.path.join(locations.representativesDir, locations.runName, locations.representativeMRNA)
        repCDSFile = os.path.join(locations.representativesDir, locations.runName, locations.representativeCDS)
        repAAFile = os.path.join(locations.representativesDir, locations.runName, locations.representativeAA)
        
        # Representatives act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "representatives",
            "-d", workDir,
            "--useGFF3", "--useBLAST", "--useSalmon"
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Representatives assert
        self.assertTrue(stderr == "", f"'representatives' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(repMRNAFile + ".ok"),
                        f"Expected '{repMRNAFile}' to have an .ok file")
        self.assertTrue(os.path.isfile(repCDSFile + ".ok"),
                        f"Expected '{repCDSFile}' to have an .ok file")
        self.assertTrue(os.path.isfile(repAAFile + ".ok"),
                        f"Expected '{repAAFile}' to have an .ok file")
        
        # DGE arrange
        rScriptFile = os.path.join(locations.dgeDir, locations.runName, locations.rScriptFile)
        salmonQCFile = os.path.join(locations.dgeDir, locations.runName, locations.salmonQCFile)
        sampleFile = os.path.join(locations.dgeDir, locations.runName, locations.salmonSampleFile)
        
        # DGE act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "dge",
            "-d", workDir
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # DGE assert
        self.assertTrue(stderr == "", f"'dge' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(rScriptFile + ".ok"),
                        f"Expected '{rScriptFile}' to have an .ok file")
        self.assertTrue(os.path.isfile(salmonQCFile + ".ok"),
                        f"Expected '{salmonQCFile}' to have an .ok file")
        self.assertTrue(os.path.isfile(sampleFile + ".ok"),
                        f"Expected '{sampleFile}' to have an .ok file")
        
        # Annotate arrange
        idMappingFile = os.path.join(dataDir, "example_id_mapping.tab")
        goOboFile = os.path.join(dataDir, "example_go.obo")
        annotFile = os.path.join(locations.annotateDir, locations.runName, locations.annotationFile)
        
        # Annotate act
        cmd = [
            PYTHON_EXE, os.path.join(baseDir, "BINge_post.py"), "annotate",
            "-d", workDir,
            "-id", idMappingFile,
            "-io", goOboFile
        ]
        returncode, stdout, stderr = run_subprocess(cmd)
        
        # Annotate assert
        self.assertTrue(stderr == "", f"'annotate' has stderr output: {stderr}")
        self.assertTrue(os.path.isfile(annotFile + ".ok"),
                        f"Expected '{annotFile}' to have an .ok file")
        
        with open(annotFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                else:
                    query, source, target_accession, gene_names, taxa_names, \
                        pct_identity, evalue, bitscore, bestwidmap, \
                        bestwgo, bestwgoplus = line.rstrip().split("\t")
                    break
        
        self.assertTrue(query == "cluster-0", f"'annotate' file should start with cluster-0, not {query}")
        self.assertTrue(target_accession.startswith("P04637"),
                        f"'annotate' file should have a best hit to P04637, not {target_accession}")       
        self.assertTrue(gene_names == "Cellular tumor antigen p53 [Cellular tumor antigen p53]",
                        f"'annotate' file should return hits to 'Cellular tumor antigen p53', not {gene_names}")
        self.assertTrue(bestwgoplus == "GO:0000001",
                        f"'annotate' file should have a gene ontology of 'GO:0000001', not {bestwgoplus}")

# Simple unit tests
class TestLocations(unittest.TestCase):
    def test_locations_init(self):
        # Arrange
        os.makedirs(workDir, exist_ok=True)
        locations = Locations(workDir)
        
        # Act
        ## N/A
        
        # Assert
        self.assertEqual(locations.workingDirectory, os.path.abspath(workDir))
    
    def test_get_gff3Files(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        expectedGff3Prefixes = set(["genome1", "genome2", "annotations1", "annotations2"])
        
        # Act
        gff3Files = locations.get_gff3Files(targetGenomes, annotatedGenomes)
        gff3Prefixes = set([ os.path.basename(x).split(".")[0] for x in gff3Files ])
        
        # Assert
        self.assertEqual(len(gff3Files), 4, "Should be 4 GFF3s after setup_working_directory()")
        self.assertEqual(expectedGff3Prefixes, gff3Prefixes, f"Should have gff3 prefixes of '{expectedGff3Prefixes}'")
    
    def test_get_sequenceFiles(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        expectedSeqPrefix = set(["genome1", "genome2", "annotations1", "annotations2", "transcriptome1"]) # no other sequences are fed in pre-extracted other than mRNA
        
        # Act
        aaFiles = locations.get_sequenceFiles(targetGenomes, annotatedGenomes, transcriptomes, "aa")
        aaPrefixes = set([ os.path.basename(x).split(".")[0] for x in aaFiles ])
        
        # Assert
        self.assertEqual(len(aaFiles), 5, "Should be 5 .aa files after setup_working_directory()")
        self.assertEqual(expectedSeqPrefix, aaPrefixes, f"Should have seq prefixes of '{expectedSeqPrefix}'")
    
    def test_resolve_runName(self):
        '''
        This test doesn't fully test the resolution process, but it's had enough real-life testing
        that this function should really just be a sanity check to make sure things didn't magically break.
        '''
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        runDir = os.path.join(locations.analysisDir, locations.runName)
        os.makedirs(runDir, exist_ok=True)
        
        # Act
        resolvedRunName = locations.resolve_runName(locations.analysisDir) # performs validation implicitly
        
        # Assert
        self.assertEqual(runDir, resolvedRunName, f"Should be '{runDir}'")

class TestTargetGenome(unittest.TestCase):
    def test_targetGenome_empty_gff3(self):
        pass ## TBD
    
    def test_targetGenome_invalid_gff3(self):
        pass ## TBD
    
    def test_targetGenome_empty_fasta(self):
        pass ## TBD
    
    def test_targetGenome_invalid_fasta(self):
        pass ## TBD
    
    def test_targetGenome_extract_1(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        with open(os.path.join(dataDir, "genome1.ttable1.aa"), "r") as fileIn:
            previousAA = fileIn.read()
        
        # Act
        targetGenome = targetGenomes[0] # just the first one
        targetGenome.extract_sequences(False, 1) # isMicrobial, translationTable
        with open(targetGenome.aa, "r") as fileIn:
            thisAA = fileIn.read()
        
        # Assert
        self.assertEqual(thisAA, previousAA, "TargetGenome extraction and genome1.ttable1.aa should be equivalent (translation table == 1)")
    
    def test_targetGenome_extract_5(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory(translationTable=5)
        with open(os.path.join(dataDir, "genome1.ttable5.aa"), "r") as fileIn:
            previousAA = fileIn.read()
        
        # Act
        targetGenome = targetGenomes[0] # just the first one
        targetGenome.extract_sequences(False, 5) # isMicrobial, translationTable
        with open(targetGenome.aa, "r") as fileIn:
            thisAA = fileIn.read()
        
        # Assert
        self.assertEqual(thisAA, previousAA, "TargetGenome extraction and genome1.ttable5.aa should be equivalent (translation table == 5)")
    
    def test_targetGenome_length_index(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        
        # Act
        targetGenome = targetGenomes[0] # just the first one
        targetGenome.length_index()
        
        # Assert
        seqLenDict = load_sequence_length_index(os.path.join(locations.genomesDir, f"{targetGenome.prefix}.lengths.pkl"))
        self.assertEqual(len(seqLenDict), 1, f"{targetGenome.prefix}.lengths.pkl should contain only 1 contig")
        self.assertEqual(seqLenDict["genome1"], 1000, "genome1 should have a length of 1000 bp")

class TestAnnotatedGenome(unittest.TestCase):
    def test_annotatedGenome_empty_gff3(self):
        pass ## TBD
    
    def test_annotatedGenome_invalid_gff3(self):
        pass ## TBD
    
    def test_annotatedGenome_empty_fasta(self):
        pass ## TBD
    
    def test_annotatedGenome_invalid_fasta(self):
        pass ## TBD
    
    def test_annotatedGenome_extract_1(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        with open(os.path.join(dataDir, "genome1.ttable1.aa"), "r") as fileIn:
            previousAA = fileIn.read()
        
        # Act
        annotatedGenome = annotatedGenomes[0] # just the first one
        annotatedGenome.extract_sequences(False, 1) # isMicrobial, translationTable
        with open(annotatedGenome.aa, "r") as fileIn:
            thisAA = fileIn.read()
        
        # Assert
        self.assertEqual(thisAA, previousAA, "AnnotatedGenome extraction and genome1.ttable1.aa should be equivalent (translation table == 1)")
    
    def test_annotatedGenome_extract_5(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory(translationTable=5)
        with open(os.path.join(dataDir, "genome1.ttable5.aa"), "r") as fileIn:
            previousAA = fileIn.read()
        
        # Act
        annotatedGenome = annotatedGenomes[0] # just the first one
        annotatedGenome.extract_sequences(False, 5) # isMicrobial, translationTable
        with open(annotatedGenome.aa, "r") as fileIn:
            thisAA = fileIn.read()
        
        # Assert
        self.assertEqual(thisAA, previousAA, "AnnotatedGenome extraction and genome1.ttable5.aa should be equivalent (translation table == 5)")

class TestTranscriptome(unittest.TestCase):
    def test_transcriptome_empty_mrna(self):
        pass ## TBD
    
    def test_transcriptome_empty_cds(self):
        pass ## TBD
    
    def test_transcriptome_empty_aa(self):
        pass ## TBD
    
    def test_transcriptome_extract_1(self):
        """
        This test wasn't what I initially expected; it isn't very useful, but it still is a behaviour that can be tested.
        Edit: However, it now raises an intentional error and hence this test no longer serves a purpose. No valid ORFs
        does not produce an empty file, it now produces an error.
        """
        pass
        
        # # Arrange
        # targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory(txomeIsPreExtracted=False)
        
        # # Act
        # transcriptome = transcriptomes[0]
        # transcriptome.extract_sequences(1) # isMicrobial
        # with open(transcriptome.aa, "r") as fileIn:
        #     thisAA = fileIn.read()
        
        # # Assert
        # self.assertEqual(thisAA, "", "Transcriptome extraction will not produce an output as no valid ORFs exist")
    
    def test_transcriptome_symlinking(self):
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        originalAA = os.path.join(dataDir, "genome3.ttable1.aa")
        
        # Act
        transcriptome = transcriptomes[0]
        transcriptome.extract_sequences(1) # isMicrobial
        aaSymlink = transcriptome.aa
        resolvedSymlink = str(Path(aaSymlink).resolve())
        
        # Assert
        self.assertTrue(os.path.exists(aaSymlink), f"'{aaSymlink}' should have been created")
        self.assertTrue(os.path.isfile(aaSymlink), f"'{aaSymlink}' should be a file")
        self.assertTrue(os.path.islink(aaSymlink), f"'{aaSymlink}' should be a symlink")
        self.assertEqual(originalAA, resolvedSymlink, f"Transcriptome.extract_sequences() should symlink back to {originalAA}")

class TestSetup(unittest.TestCase):
    def test_json_fidelity_1(self):
        '''
        Ensure that inputs_to_json() and json_to_inputs() maintain argument fidelity
        '''
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        isMicrobial = True
        
        targetKeys = targetGenomes[0].__dict__.keys()
        annotatedKeys = annotatedGenomes[0].__dict__.keys()
        txomeKeys = transcriptomes[0].__dict__.keys()
        
        # Act
        inputs_to_json(locations, targetGenomes, annotatedGenomes, transcriptomes, isMicrobial)
        _targetGenomes, _annotatedGenomes, _transcriptomes, _isMicrobial = json_to_inputs(locations)
        
        targetsSame = [ targetGenomes[i].__dict__[key] == _targetGenomes[i].__dict__[key] for i in range(len(targetGenomes)) for key in targetKeys ]
        annotatedSame = [ annotatedGenomes[i].__dict__[key] == _annotatedGenomes[i].__dict__[key] for i in range(len(annotatedGenomes)) for key in annotatedKeys ]
        txomeSame = [ transcriptomes[i].__dict__[key] == _transcriptomes[i].__dict__[key] for i in range(len(transcriptomes)) for key in txomeKeys ]
        
        # Assert
        self.assertTrue(all(targetsSame), "JSON reconstitution of TargetGenomes should give identical values")
        self.assertTrue(all(annotatedSame), "JSON reconstitution of AnnotatedGenomes should give identical values")
        self.assertTrue(all(txomeSame), "JSON reconstitution of Transcriptomes should give identical values")
        self.assertEqual(isMicrobial, _isMicrobial, f"JSON reconstitution of isMicrobial should give {isMicrobial}")
    
    def test_json_fidelity_2(self):
        '''
        Ensure that inputs_to_json() and json_to_inputs() maintain argument fidelity
        '''
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory(
            targetHasGFF3=False,
            txomeIsPreExtracted=True # this was originally False, but it now triggers an (intentional) error
        )
        isMicrobial = False
        
        targetKeys = targetGenomes[0].__dict__.keys()
        annotatedKeys = annotatedGenomes[0].__dict__.keys()
        txomeKeys = transcriptomes[0].__dict__.keys()
        
        # Act
        inputs_to_json(locations, targetGenomes, annotatedGenomes, transcriptomes, isMicrobial)
        _targetGenomes, _annotatedGenomes, _transcriptomes, _isMicrobial = json_to_inputs(locations)
        
        targetsSame = [ targetGenomes[i].__dict__[key] == _targetGenomes[i].__dict__[key] for i in range(len(targetGenomes)) for key in targetKeys ]
        annotatedSame = [ annotatedGenomes[i].__dict__[key] == _annotatedGenomes[i].__dict__[key] for i in range(len(annotatedGenomes)) for key in annotatedKeys ]
        txomeSame = [ transcriptomes[i].__dict__[key] == _transcriptomes[i].__dict__[key] for i in range(len(transcriptomes)) for key in txomeKeys ]
        
        # Assert
        self.assertTrue(all(targetsSame), "JSON reconstitution of TargetGenomes should give identical values")
        self.assertTrue(all(annotatedSame), "JSON reconstitution of AnnotatedGenomes should give identical values")
        self.assertTrue(all(txomeSame), "JSON reconstitution of Transcriptomes should give identical values")
        self.assertEqual(isMicrobial, _isMicrobial, f"JSON reconstitution of isMicrobial should give {isMicrobial}")

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

class TestValidation(unittest.TestCase):
    def test_check_for_duplicates_1(self):
        'This check should trigger an exception as the target and annotated genomes have duplicated IDs'
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        
        # Act & Assert
        with self.assertRaises(ValueError):
            check_for_duplicates(locations.get_sequenceFiles(
                targetGenomes, annotatedGenomes,transcriptomes, "aa")
            )
    
    def test_check_for_duplicates_2(self):
        'This check should NOT trigger an exception as there are no duplicated IDs'
        # Arrange
        file1 = os.path.join(dataDir, "genome1.ttable1.aa")
        file2 = os.path.join(dataDir, "genome3.ttable1.aa")
        
        # Act & Assert
        try:
            check_for_duplicates([file1, file2])
            self.assertTrue(True)
        except:
            self.assertTrue(False, f"Duplication issue should not exist for '{file1}' and '{file2}'")
    
    def test_check_for_consistency_1(self):
        'This check should trigger an exception as the provided files have different sequence identifiers'
        # Arrange
        mrna = os.path.join(dataDir, "genome1.ttable1.mrna")
        cds = os.path.join(dataDir, "genome1.ttable1.cds")
        aa = os.path.join(dataDir, "genome3.ttable1.aa")
        
        # Act & Assert
        with self.assertRaises(ValueError):
            check_for_seqid_consistency(mrna, cds, aa)
    
    def test_check_for_consistency_2(self):
        'This check should NOT trigger an exception as the provided files all have consistent sequence identifiers'
        # Arrange
        targetGenomes, annotatedGenomes, transcriptomes, locations = setup_working_directory()
        
        # Act
        transcriptome = transcriptomes[0]
        mrna, cds, aa = transcriptome.mrna, transcriptome.cds, transcriptome.aa
        
        # Act & Assert
        try:
            check_for_seqid_consistency(mrna, cds, aa)
            self.assertTrue(True)
        except:
            self.assertTrue(False, f"Consistency issue should not exist for '{transcriptome}'")

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

class TestGff3ToFasta(unittest.TestCase):
    def test_ncrna_gff3_1(self):
        'No sequences should be created, and we consequently expect an error to occur within gff3_to_fasta()'
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        make_working_directory(locations)
        
        gff3FileIn = os.path.join(dataDir, "no_gene_features.gff3")
        fastaFileIn = os.path.join(dataDir, "genome1.fasta")
        
        mrnaFileOut = os.path.join(locations.sequencesDir, "test.mrna")
        cdsFileOut = os.path.join(locations.sequencesDir, "test.cds")
        protFileOut = os.path.join(locations.sequencesDir, "test.aa")
        
        isMicrobial = False
        translationTable = 1
        
        # Act & Assert
        with self.assertRaises(ValueError):
            gff3_to_fasta(gff3FileIn, fastaFileIn, mrnaFileOut, cdsFileOut,
                          protFileOut, isMicrobial, translationTable)
    
    def test_ncrna_gff3_2(self):
        'Similar to test_ncrna_gff3_1(), but with isMicrobial==True'
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        make_working_directory(locations)
        
        gff3FileIn = os.path.join(dataDir, "no_gene_features.gff3")
        fastaFileIn = os.path.join(dataDir, "genome1.fasta")
        
        mrnaFileOut = os.path.join(locations.sequencesDir, "test.mrna")
        cdsFileOut = os.path.join(locations.sequencesDir, "test.cds")
        protFileOut = os.path.join(locations.sequencesDir, "test.aa")
        
        isMicrobial = True
        translationTable = 11
        
        # Act & Assert
        with self.assertRaises(ValueError):
            gff3_to_fasta(gff3FileIn, fastaFileIn, mrnaFileOut, cdsFileOut,
                          protFileOut, isMicrobial, translationTable)

class TestTxomeToOrfs(unittest.TestCase):
    def test_txome_to_orfs_1(self):
        'No sequences should be created, and we consequently expect an error to occur within txome_to_orfs()'
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        make_working_directory(locations)
        
        mrnaFileIn = os.path.join(dataDir, "genome3.ttable1.mrna")
        
        cdsFileOut = os.path.join(locations.sequencesDir, "test.cds")
        protFileOut = os.path.join(locations.sequencesDir, "test.aa")
        
        translationTable = 1
        
        # Act & Assert
        with self.assertRaises(ValueError):
            txome_to_orfs(mrnaFileIn, cdsFileOut, protFileOut, translationTable)
    
    def test_txome_to_orfs_2(self):
        'Similar to test_txome_to_orfs_1(), but with translationTable==11'
        # Arrange: obtain variables
        cleanup_working_directory()
        locations = Locations(workDir)
        make_working_directory(locations)
        
        mrnaFileIn = os.path.join(dataDir, "genome3.ttable1.mrna")
        
        cdsFileOut = os.path.join(locations.sequencesDir, "test.cds")
        protFileOut = os.path.join(locations.sequencesDir, "test.aa")
        
        translationTable = 11
        
        # Act & Assert
        with self.assertRaises(ValueError):
            txome_to_orfs(mrnaFileIn, cdsFileOut, protFileOut, translationTable)

class TestGmapGFF3(unittest.TestCase):
    def test_iterate_through_empty_file(self):
        'GmapGFF3 class should error out on an empty file'
        # Arrange
        gmapFile = os.path.join(dataDir, "empty_file")
        gmapGff3Obj = GmapGFF3(gmapFile)
        
        # Act and assert
        with self.assertRaises(ValueError):
            for feature in gmapGff3Obj:
                pass
    
    def test_iterate_through_gff3(self):
        'GmapGFF3 should parse out the expected information from a normal GMAP file'
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
        self.assertEqual(numFeatures, 2, "GMAP GFF3 should contain 2 features")
        self.assertEqual(mrnaIDs, ['genome1.1', 'genome1.2'],
                         "Should be ['genome1.1', 'genome1.2']")

class TestGFF3Graph(unittest.TestCase):
    def test_init_with_empty_file(self):
        'GFF3Graph class should error out during initial parsing of an empty file'
        # Arrange
        gff3File = os.path.join(dataDir, "empty_file")
        
        # Act and assert
        with self.assertRaises(ValueError):
            graphObj = GFF3Graph(gff3File)
    
    def test_iterate_through_gff3(self):
        'GFF3Graph should parse out the expected information from a normal GFF3 file'
        # Arrange
        gff3File = os.path.join(dataDir, "gmap_normal.gff3")
        graphObj = GFF3Graph(gff3File)
        
        # Act
        numFeatures = 0
        mrnaIDs = []
        for featureID in graphObj.ftypes["mRNA"]:
            feature = graphObj[featureID]
            numFeatures += 1
            mrnaIDs.append(feature.ID)
        
        # Assert
        self.assertEqual(numFeatures, 2, "GMAP GFF3 should contain 2 features")
        self.assertEqual(mrnaIDs, ['genome1.1.mrna1', 'genome1.2.mrna1'],
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
