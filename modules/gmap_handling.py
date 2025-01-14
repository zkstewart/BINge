import os, sys, re
from .fasta_handling import remove_sequence_from_fasta
from .thread_workers import BasicProcess
from .validation import touch_ok

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_MapIO import GMAP_DB, GMAP

# Multithreaded functions and classes
class GmapIndexProcess(BasicProcess):
    '''
    Handles GMAP indexing in a separate thread.
    
    Parameters:
        fasta -- a string indicating the location of a FASTA file for GMAP indexing.
        gmapDir -- a string indicating the location of the GMAP executable files.
    '''
    def task(self, fasta, gmapDir):
        db = GMAP_DB(fasta, gmapDir)
        if not db.index_exists():
            db.index()

# Other functions
def setup_gmap_indices(locations, gmapDir, threads):
    '''
    Will take the genome files in the workingDirectory 'genomes' subdir and
    generate indices for all files.
    
    Parameters:
        locations -- a Locations object with attributes for directory locations.
        gmapDir -- a string indicating the location where 'gmap' and 'gmap_build' are found.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Locate subdirectory containing files
    if not os.path.isdir(locations.genomesDir):
        raise FileNotFoundError(f"setup_gmap_indices failed because '{locations.genomesDir}' doesn't exist somehow?")
    
    # Locate all genome files for indexing
    genomeFiles = [
        [ os.path.join(locations.genomesDir, f), f.split(".fasta")[0] ]
        for f in os.listdir(locations.genomesDir)
        if f.endswith(".fasta")
    ]
    if not len(genomeFiles) > 0:
        raise FileNotFoundError(f"setup_gmap_indices failed because '{locations.genomesDir}' is empty somehow?")
    
    # Narrow down files to those needing indexing
    needsIndexing = []
    for genomeFile, genomePrefix in genomeFiles:
        db = GMAP_DB(genomeFile, gmapDir)
        if not db.index_exists():
            needsIndexing.append(genomeFile)
    
    # Process each db in need of indexing via threading
    if len(needsIndexing) > 0:
        print(f"# Setting up GMAP indices...")
        
        for i in range(0, len(needsIndexing), threads): # only process n (threads) files at a time
            processing = []
            for x in range(threads): # begin processing n files
                if i+x < len(needsIndexing): # parent loop may excess if n > the number of files needing indexing
                    genomeFile = needsIndexing[i+x]
                    
                    indexWorkerThread = GmapIndexProcess(genomeFile, gmapDir)
                    indexWorkerThread.start()
                    processing.append(indexWorkerThread)
            
            # Gather results
            for indexWorkerThread in processing:
                indexWorkerThread.join()
                indexWorkerThread.check_errors()

def auto_gmapping(locations, gmapDir, threads):
    '''
    Will take the transcriptome and annotation files in workingDirectory and
    use GMAP to align them to every file in the 'genomes' subdirectory.
    
    Parameters:
        locations -- a Locations object with attributes for directory locations.
        gmapDir -- a string indicating the location where 'gmap' and 'gmap_build' are found.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    PROBLEM_REGEX = re.compile(r"Problem sequence: (.+?) \(.+?\)\n")
    
    # Create subdirectory for output files (if not already existing)
    os.makedirs(locations.mappingDir, exist_ok=True)
    
    # Locate subdirectory containing files
    if not os.path.isdir(locations.genomesDir):
        raise FileNotFoundError(f"auto_gmapping failed because '{locations.genomesDir}' doesn't exist somehow?")
    
    # Locate all genome files for alignment target
    genomeFiles = [
        [ os.path.join(locations.genomesDir, f), f.split(".fasta", maxsplit=1)[0] ]
        for f in os.listdir(locations.genomesDir)
        if f.endswith(".fasta")
    ]
    if not len(genomeFiles) > 0:
        raise FileNotFoundError(f"auto_gmapping failed because '{locations.genomesDir}' is empty somehow?")
    
    # Locate all sequence files for alignment query
    queryFiles = [
        [ os.path.join(locations.sequencesDir, f), f.split(".cds", maxsplit=1)[0] ]
        for f in os.listdir(locations.sequencesDir)
        if f.endswith(".cds")
    ] + [
        [ os.path.join(locations.genomesDir, f), f.split(".cds", maxsplit=1)[0] ]
        for f in os.listdir(locations.genomesDir)
        if f.endswith(".cds")
    ]
    
    if not len(queryFiles) > 0:
       raise FileNotFoundError(f"auto_gmapping failed because '{locations.sequencesDir}' contains no query files somehow?")
    
    # Iteratively perform GMAP search for all combinations
    notifiedOnce = False
    for queryFile, queryPrefix in queryFiles:
        originalQuery = queryFile # remember what the original query file was if we end up using tmp files
        tmpQueryFile = f"{originalQuery}.tmp" # temporary file for problem sequences
        
        for genomeFile, genomePrefix in genomeFiles:
            outputFileName = os.path.join(locations.mappingDir, f"{queryPrefix}_to_{genomePrefix}_gmap.gff3")
            if not os.path.exists(outputFileName) or not os.path.exists(outputFileName + ".ok"):
                # Give user a heads up
                if not notifiedOnce:
                    print(f"# Running GMAP alignment...")
                    notifiedOnce = True
                
                # Continue processing
                problemIDs = []
                while True:
                    try:
                        # Use a temporary file if we have problem sequences
                        if os.path.exists(tmpQueryFile):
                            queryFile = tmpQueryFile
                        else:
                            queryFile = originalQuery
                        
                        # Run GMAP
                        gmapper = GMAP(queryFile, genomeFile, gmapDir, threads)
                        assert gmapper.index_exists(), \
                            f"auto_gmapping failed because '{genomeFile}' doesn't have an index?"
                        
                        gmapper.gmap(outputFileName, force=True) # allows overwrite if .ok file is missing
                        
                        # Clean up temporary file (if it exists) on successful run
                        if os.path.exists(f"{originalQuery}.tmp"):
                            os.unlink(f"{originalQuery}.tmp")
                    
                    except Exception as e:
                        # Delete the failed file
                        if os.path.exists(outputFileName):
                            os.unlink(outputFileName)
                        
                        # See if this error is because of a problem sequence
                        if "Problem sequence" in e.args[0]:
                            # Identify the sequence ID
                            problemIDs.append(PROBLEM_REGEX.search(e.args[0]).groups()[0])
                            
                            # Remove the sequence from the file
                            remove_sequence_from_fasta(originalQuery, problemIDs, tmpQueryFile, force=True)
                            
                            # Try again
                            continue
                        else:
                            raise e # re-raise the error if it's not a problem sequence error
                    
                    # If no errors, then exit out of the while loop
                    break
            touch_ok(outputFileName)
