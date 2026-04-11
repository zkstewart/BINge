import os, sys, re, subprocess, platform
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from fasta_handling import remove_sequence_from_fasta
from thread_workers import BasicProcess
from validation import touch_ok

class DirectoryNotFoundError(Exception):
    pass

class GMAP:
    '''
    The GMAP Class provides easy access to the GMAP search function to perform
    mapping of FASTAs against each other. Object values can be set to control
    the parameters specified for a GMAP search.
    
    Attributes:
        query -- a string indicating the location of a FASTA file.
        db -- a string indicating the location of a GMAP index directory.
        gmapDir -- a string indicating the location of GMAP binaries.
        threads -- (OPTIONAL) an integer value specifying the number of threads
                   to run when performing GMAP search. Defaults to 1.
    '''
    def __init__(self, query, db, gmapDir, threads=1):
        self.query = query
        self.gmapDir = gmapDir
        self.db = db
        self.threads = threads
        
        self.isGMAP = True # flag to check object type
        
        # Set default values
        self.outputFormat = 2
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
            raise FileNotFoundError("GMAP query parameter does not point to an existing file")
        
        self._query = os.path.abspath(value)
    
    @property
    def db(self):
        return self._db
    
    @db.setter
    def db(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isdir(value):
            raise DirectoryNotFoundError("GMAP db parameter does not point to an existing directory")
        
        self._db = os.path.abspath(value)
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise DirectoryNotFoundError("gmapDir does not point to an existing directory")
        
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
            raise FileNotFoundError("gmapExe does not point to an existing file")
        
        self._gmapExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise ValueError("Number of threads must be more than 0")
        
        self._threads = value
    
    @property
    def outputFormat(self):
        return self._outputFormat
    
    @outputFormat.setter
    def outputFormat(self, value):
        assert isinstance(value, int)
        if value < 1 or value > 9:
            raise ValueError("GMAP --format only supports output formats 1-9")
        
        self._outputFormat = value
    
    @property
    def npaths(self):
        return self._npaths
    
    @npaths.setter
    def npaths(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise ValueError("GMAP --npaths must be >= 0")
        
        self._npaths = value
    
    @property
    def chimeraMargin(self):
        return self._chimeraMargin
    
    @chimeraMargin.setter
    def chimeraMargin(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise ValueError("GMAP --chimera-margin must be >= 0")
        
        self._chimeraMargin = value
    
    @property
    def batch(self):
        return self._batch
    
    @batch.setter
    def batch(self, value):
        assert isinstance(value, int)
        if value < 0 or value > 5:
            raise ValueError("GMAP --batch must be in the range 0-5")
        
        self._batch = value
    
    @property
    def maxIntronLengthMiddle(self):
        return self._maxIntronLengthMiddle
    
    @maxIntronLengthMiddle.setter
    def maxIntronLengthMiddle(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise ValueError("GMAP --max-intronlength-middle must be a positive integer")
        
        self._maxIntronLengthMiddle = value
    
    @property
    def maxIntronLengthEnds(self):
        return self._maxIntronLengthEnds
    
    @maxIntronLengthEnds.setter
    def maxIntronLengthEnds(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise ValueError("GMAP --max-intronlength-ends must be a positive integer")
        
        self._maxIntronLengthEnds = value
    
    def gmap(self, outFile):
        '''
        Performs the GMAP search operation.
        
        Parameters:
            outFile -- a string indicating the location to write the output file to. File must not already exist!
        '''
        # Format command
        if self.query.endswith(".gz"):
            cmd = ["zcat", self.query, "|"]
        else:
            cmd = ["cat", self.query, "|"]
        cmd.extend([
            self.gmapExe, "-D", os.path.dirname(self.db), "-d", os.path.basename(self.db),
            "-f", str(self.outputFormat), "-n", str(self.npaths), "-x", str(self.chimeraMargin),
            "-B", str(self.batch), "-t", str(self.threads),
            f"--max-intronlength-middle={self.maxIntronLengthMiddle}",
            f"--max-intronlength-ends={self.maxIntronLengthEnds}",
            ">", outFile
        ])
        
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

def build_index(buildExe, fastaFile, dbDirectory):
    '''
    Parameters:
        buildExe -- a string indicating the full path to the 'gmap_build' executable.
        fastaFile -- a string pointing to a FASTA file to be indexed.
        dbDirectory -- a string pointing to the full location to write the index.
    '''
    # Format command
    cmd = [
        buildExe, "-D", os.path.dirname(dbDirectory), # -D chooses the directory to write the index
        "-d", f"{os.path.basename(dbDirectory)}" # -d indicate the name for the index output
    ]
    if fastaFile.endswith(".gz"):
        cmd.append("-g") # -g tells GMAP to gunzip the file
    cmd.append(fastaFile)
    
    # Run indexing
    if platform.system() != "Windows":
        run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        run_index = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    indexout, indexerr = run_index.communicate()
    
    # Raise exception if final part of stderr from successful run isn't found
    indexerr = indexerr.decode("utf-8")
    if (not "Writing localdb sarrays" in indexerr) and (not indexerr.rstrip("\r\n ").endswith("Done")):
        raise Exception('GMAP indexing error text below\n' + indexerr)

def auto_gmapping(targetGenomes, annotatedGenomes, transcriptomes, mappingDir, gmapDir, threads):
    '''
    Will take the transcriptome and annotation files in workingDirectory and
    use GMAP to align them to every file in the 'genomes' subdirectory.
    
    This function may potentially go through an elaborate process using temporary query files
    and remove_sequence_from_fasta(). This can be invoked where a GMAP error occurs with the program
    noting which sequence caused the program to crash. More recent GMAP versions appear to have
    few errors like this, but it used to be a common occurrence in the early phase of BINge's development.
    It remains in place on the off chance it helps to resolve an issue.
    
    Parameters:
        targetGenomes -- a list of TargetGenome objects, to be used as references
        annotatedGenomes -- a list of AnnotatedGenome objects, to be queried to reference;
                            this can be empty
        transcriptomes -- a list of Transcriptome objects, to be queried to reference;
                          this can be empty
        mappingDir -- a string pointing to the location to write GMAP results to.
        gmapDir -- a string indicating the location where 'gmap' and 'gmap_build' are found.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    PROBLEM_REGEX = re.compile(r"Problem sequence: (.+?) \(.+?\)\n")
    
    # Locate all sequence files for alignment query
    queries = [*[tg for tg in targetGenomes if tg.cds != None], # .cds is set if target/ref genome was provided with a GFF3
               *annotatedGenomes,
               *transcriptomes]
    
    # Iteratively perform GMAP search for all combinations
    for query in queries:        
        originalQuery = query.cds # remember what the original query file was if we end up using tmp files
        tmpQueryFile = f"{originalQuery}.tmp" # temporary file for problem sequences
        
        for genome in targetGenomes:
            outputFileName = os.path.join(mappingDir, f"{query.prefix}_to_{genome.prefix}_gmap.gff3")
            if not (os.path.exists(outputFileName) and os.path.exists(outputFileName + ".ok")):
                print(f"# Performing GMAP alignment of '{query.prefix}' to '{genome.prefix}'")
                
                problemIDs = []
                while True:
                    try:
                        # Use a temporary file if we have problem sequences
                        if os.path.exists(tmpQueryFile):
                            queryFile = tmpQueryFile
                        else:
                            queryFile = originalQuery
                        
                        # Run GMAP
                        #db = genomeFile.split(".fasta")[0] + ".gmap" # temporary hack, do not leave this in forever!
                        #db = os.path.join(genome.directory, f"{genome.prefix}.gmap")
                        gmapper = GMAP(queryFile, genome.gmap_db, gmapDir, threads)
                        gmapper.gmap(outputFileName)
                        
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
