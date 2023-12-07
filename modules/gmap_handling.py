import os, sys

from .gff3_handling import GFF3
from .thread_workers import GmapIndexThread

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_MapIO import GMAP_DB, GMAP

def setup_gmap_indices(workingDirectory, gmapDir, threads):
    '''
    Will take the genome files in the workingDirectory 'genomes' subdir and
    generate indices for all files.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
        gmapDir -- a string indicating the location where 'gmap' and 'gmap_build' are found.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Locate subdirectory containing files
    genomesDir = os.path.join(workingDirectory, "genomes")
    assert os.path.isdir(genomesDir), \
        f"setup_gmap_indices failed because '{genomesDir}' doesn't exist somehow?"
    
    # Locate all genome files for indexing
    genomeFiles = [
        [ os.path.join(genomesDir, f), f.split(".fasta")[0] ]
        for f in os.listdir(genomesDir)
        if f.endswith(".fasta")
    ]
    assert len(genomeFiles) > 0, \
        f"setup_gmap_indices failed because '{genomesDir}' is empty somehow?"
    
    # Narrow down files to those needing indexing
    needsIndexing = []
    for genomeFile, genomePrefix in genomeFiles:
        db = GMAP_DB(genomeFile, gmapDir)
        if not db.index_exists():
            needsIndexing.append(genomeFile)
    
    # Process each db in need of indexing via threading
    for i in range(0, len(needsIndexing), threads): # only process n (threads) files at a time
        processing = []
        for x in range(threads): # begin processing n files
            if i+x < len(needsIndexing): # parent loop may excess if n > the number of files needing indexing
                genomeFile = needsIndexing[i+x]
                
                indexWorkerThread = GmapIndexThread(genomeFile, gmapDir)
                processing.append(indexWorkerThread)
                indexWorkerThread.start()
        
        # Gather results
        for indexWorkerThread in processing:
            indexWorkerThread.join()

def auto_gmapping(workingDirectory, gmapDir, threads):
    '''
    Will take the transcriptome and annotation files in workingDirectory and
    use GMAP to align them to every file in the 'genomes' subdirectory.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
        gmapDir -- a string indicating the location where 'gmap' and 'gmap_build' are found.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Create subdirectory for output files (if not already existing)
    mappingDir = os.path.join(workingDirectory, "mapping")
    os.makedirs(mappingDir, exist_ok=True)
    
    # Locate subdirectory containing files
    genomesDir = os.path.join(workingDirectory, "genomes")
    assert os.path.isdir(genomesDir), \
        f"auto_gmapping failed because '{genomesDir}' doesn't exist somehow?"
    
    # Locate all genome files for alignment target
    genomeFiles = [
        [ os.path.join(genomesDir, f), f.split(".fasta")[0] ]
        for f in os.listdir(genomesDir)
        if f.endswith(".fasta")
    ]
    assert len(genomeFiles) > 0, \
        f"auto_gmapping failed because '{genomesDir}' is empty somehow?"
    
    # Locate all sequence files for alignment query
    queryFiles = [
        [ os.path.join(workingDirectory, f), f.split(".nucl")[0] ]
        for f in os.listdir(workingDirectory)
        if f.endswith(".nucl")
    ]
    
    # Iteratively perform GMAP search for all combinations
    gmapFiles = []
    for queryFile, queryPrefix in queryFiles:
        for genomeFile, genomePrefix in genomeFiles:
            outputFileName = os.path.join(mappingDir, f"{queryPrefix}_to_{genomePrefix}_gmap.gff3")
            if not os.path.exists(outputFileName):
                gmapper = GMAP(queryFile, genomeFile, gmapDir, threads)
                assert gmapper.index_exists(), \
                    f"auto_gmapping failed because '{genomeFile}' doesn't have an index?"
                
                gmapper.gmap(outputFileName)
            gmapFiles.append(outputFileName)
    
    return gmapFiles
