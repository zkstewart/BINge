import os, sys, re, subprocess, platform
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from fasta_handling import remove_sequence_from_fasta
from thread_workers import BasicProcess
from validation import touch_ok

class GMAP_DB:
    '''
    The GMAP_DB Class encapsulates the logic of indexing a FASTA file using GMAP.
    
    Attributes:
        fasta (REQUIRED) -- a string indicating the location of a FASTA file.
        gmapDir (REQUIRED) -- a string indicating the location of GMAP binaries.
    '''
    def __init__(self, fasta, gmapDir):
        self.fasta = fasta
        self.gmapDir = gmapDir
        
        self.isGMAP_DB = True # flag to check object type
    
    @property
    def fasta(self):
        return self._fasta
    
    @fasta.setter
    def fasta(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Fasta parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._fasta = value
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("gmapDir does not point to an existing directory"))
        
        self._gmapDir = value
        self.buildExe = os.path.join(value, "gmap_build")
    
    @property
    def buildExe(self):
        return self._buildExe
    
    @buildExe.setter
    def buildExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("buildExe does not point to an existing file"))
        
        self._buildExe = value
    
    def index_exists(self):
        '''
        Relies on a simple assumption that a .gmap directory's presence
        indicates that a GMAP database was successfully created
        from the FASTA file.
        
        Returns:
            dbExists -- a Boolean where True means the database exists,
                        and False means it does not exist.
        '''
        return os.path.isdir(f"{self.fasta}.gmap")
    
    def index(self):
        '''
        Makes a database out of the .fasta value for use in BLAST search.
        '''
        # Skip if index exists
        if self.index_exists():
            raise FileExistsError("GMAP index already exists!")
        
        # Specify file locations
        fasta = os.path.abspath(self.fasta)
        
        # Format command
        cmd = [
            self.buildExe, "-D", os.path.dirname(fasta),
            "-d", f"{os.path.basename(fasta)}.gmap", fasta
        ]
        
        # Run indexing
        if platform.system() != "Windows":
            run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_index = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        
        # Raise exception if final part of stderr from successful run isn't found
        if "Writing localdb sarrays" not in indexerr.decode("utf-8"):
            raise Exception('GMAP indexing error text below\n' + indexerr.decode("utf-8"))

class GMAP:
    '''
    The GMAP Class provides easy access to the GMAP search function to perform
    mapping of FASTAs against each other. Object values can be set to control
    the parameters specified for a GMAP search.
    
    Attributes:
        query (REQUIRED) -- a string indicating the location of a FASTA file.
        target (REQUIRED) -- a string indicating the location of a FASTA file.
        gmapDir (REQUIRED) -- a string indicating the location of GMAP binaries.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing BLAST search. Defaults to 1.
    '''
    def __init__(self, query, target, gmapDir, threads=1):
        self.query = query
        self.gmapDir = gmapDir
        self.target = target # sets self.db
        self.threads = threads
        
        self.isGMAP = True # flag to check object type
        
        # Set default values
        self.outputFormat = 2 # my default, not GMAPs
        self.npaths = 5
        self.chimeraMargin = 30
        self.batch = 2
        self.maxIntronLengthMiddle = 500000
        self.maxIntronLengthEnds = 10000
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Query parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._query = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Target parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._target = value
        self.db = GMAP_DB(value, self.gmapDir)
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("gmapDir does not point to an existing directory"))
        
        self._gmapDir = value
        self.gmapExe = os.path.join(value, "gmap")
    
    @property
    def gmapExe(self):
        return self._gmapExe
    
    @gmapExe.setter
    def gmapExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("gmapExe does not point to an existing file"))
        
        self._gmapExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("Number of threads must be more than 0")
        
        self._threads = value
    
    @property
    def outputFormat(self):
        return self._outputFormat
    
    @outputFormat.setter
    def outputFormat(self, value):
        assert isinstance(value, int)
        if value < 1 or value > 9:
            raise Exception("GMAP --format only supports output formats 1-9")
        
        self._outputFormat = value
    
    @property
    def npaths(self):
        return self._npaths
    
    @npaths.setter
    def npaths(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("GMAP --npaths must be >= 0")
        
        self._npaths = value
    
    @property
    def chimeraMargin(self):
        return self._chimeraMargin
    
    @chimeraMargin.setter
    def chimeraMargin(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("GMAP --chimera-margin must be >= 0")
        
        self._chimeraMargin = value
    
    @property
    def batch(self):
        return self._batch
    
    @batch.setter
    def batch(self, value):
        assert isinstance(value, int)
        if value < 0 or value > 5:
            raise Exception("GMAP --batch must be in the range 0-5")
        
        self._batch = value
    
    @property
    def maxIntronLengthMiddle(self):
        return self._maxIntronLengthMiddle
    
    @maxIntronLengthMiddle.setter
    def maxIntronLengthMiddle(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise Exception("GMAP --max-intronlength-middle must be a positive integer")
        
        self._maxIntronLengthMiddle = value
    
    @property
    def maxIntronLengthEnds(self):
        return self._maxIntronLengthEnds
    
    @maxIntronLengthEnds.setter
    def maxIntronLengthEnds(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise Exception("GMAP --max-intronlength-ends must be a positive integer")
        
        self._maxIntronLengthEnds = value
    
    def index_exists(self):
        '''
        Checks if the .target value is in a GMAP index/db yet.
        '''
        return self.db.index_exists()
    
    def index(self):
        '''
        Makes a database out of the .target value for use in BLAST search.
        '''
        # Skip if index exists
        if self.db.index_exists():
            raise FileExistsError("GMAP index already exists!")
        
        # Create index otherwise
        self.db.index()
    
    def gmap(self, outFile, force=False):
        '''
        Performs the GMAP search operation.
        
        This method does not use .query and .target because they may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            outFile -- a string indicating the location to write the output file to. File must not already exist!
            force -- a boolean indicating whether to overwrite the output file if it already exists.
        '''
        # Skip if index does not exist
        if not self.db.index_exists():
            raise FileNotFoundError("GMAP index does not exist! Cannot perform search.")
        
        # Error if output file already exists
        if os.path.isfile(outFile) and force == False:
            raise FileExistsError(f"GMAP output file '{outFile}' already exists! Will not overwrite.")
        
        # Specify file locations
        query = os.path.abspath(self.query)
        target = os.path.abspath(self.target)
        
        # Format command
        cmd = [
            self.gmapExe, "-D", os.path.dirname(target), "-d", f"{os.path.basename(target)}.gmap",
            "-f", str(self.outputFormat), "-n", str(self.npaths), "-x", str(self.chimeraMargin),
            "-B", str(self.batch), "-t", str(self.threads),
            f"--max-intronlength-middle={self.maxIntronLengthMiddle}",
            f"--max-intronlength-ends={self.maxIntronLengthEnds}",
            query, ">", outFile
        ]
        
        # Perform GMAP search
        if platform.system() != "Windows":
            run_gmap = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_gmap = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        gmapout, gmaperr = run_gmap.communicate()
        if "queries/sec)" not in gmaperr.decode("utf-8"):
            raise Exception('GMAP searching error text below\n' +
                            gmaperr.decode("utf-8"))

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
                try:
                    indexWorkerThread.check_errors()
                except Exception as e:
                    if not str(e).rstrip("\r\n ").endswith("Done"):
                        raise e

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
    queryFiles = locations.get_sequenceFiles(".cds")
    
    # Iteratively perform GMAP search for all combinations
    notifiedOnce = False
    for queryFile in queryFiles:
        queryPrefix = os.path.basename(queryFile).split(".cds", maxsplit=1)[0]
        
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
