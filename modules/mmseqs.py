import os, sys, subprocess
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from validation import touch_ok

class MM_DB:
    '''
    The MM_DB Class provides the logic for creating and handling MMSeqs2 databases.
    
    Attributes:
        fasta (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        mmseqsDir (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                is found.
        tmpDir (REQUIRED) -- a string location for where MMSeqs2 should keep temp files.
        molecule (REQUIRED) -- a string indicating the type of molecule the database is for;
                               acceptable values are 'protein' or 'nucleotide'.
        threads (OPTIONAL) -- a positive integer for how many threads to use when indexing
                              a database file.
    '''
    def __init__(self, fasta, mmseqsDir, tmpDir, molecule, threads=1):
        self.fasta = fasta
        self.mmseqsDir = mmseqsDir
        self.tmpDir = tmpDir
        self.molecule = molecule
        self.threads = threads
        
        self.isMM_DB = True # flag to check object type
    
    @property
    def fasta(self):
        return self._fasta
    
    @fasta.setter
    def fasta(self, value):
        assert type(value).__name__ == "str"
        if type(value).__name__ == "str" and not os.path.isfile(value):
            raise Exception(("Fasta parameter does not point to an existing file location"))
        
        self._fasta = os.path.abspath(value)
    
    @property
    def mmseqsDir(self):
        return self._mmseqsDir
    
    @mmseqsDir.setter
    def mmseqsDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("mmseqsDir does not point to an existing directory"))
        
        self._mmseqsDir = os.path.abspath(value)
        self.mmseqsExe = os.path.join(self._mmseqsDir, "mmseqs")
    
    @property
    def mmseqsExe(self):
        return self._mmseqsExe
    
    @mmseqsExe.setter
    def mmseqsExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("mmseqsExe does not point to an existing file"))
        
        self._mmseqsExe = value
    
    @property
    def tmpDir(self):
        return self._tmpDir
    
    @tmpDir.setter
    def tmpDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        value = os.path.abspath(value)
        
        if not os.path.isdir(os.path.dirname(value)):
            raise Exception((f"tmpDir's parent location ('{os.path.dirname(value)}') " +
                             "does not exist"))
        
        if not os.path.isdir(value):
            os.mkdir(value)
        
        self._tmpDir = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    @property
    def molecule(self):
        return self._molecule
    
    @molecule.setter
    def molecule(self, value):
        assert isinstance(value, str)
        assert value in ["protein", "nucleotide"], \
            "molecule must be either 'protein' or 'nucleotide'"
        
        self._molecule = value
    
    @property
    def db(self):
        return f"{self.fasta}_seqDB" # self.fasta is already abspath
    
    @property
    def db_basename(self):
        return os.path.basename(self.db)
    
    @property
    def db_dirname(self):
        return os.path.dirname(self.db)
    
    @property
    def db_ok(self):
        return self.db + ".db.ok"
    
    @property
    def index_ok(self):
        return self.db + ".index.ok"
    
    def generate(self):
        '''
        Creates a sequence DB for MMSeqs2 if it does not already exist.
        '''
        # Skip if DB already exists
        if os.path.isfile(self.db) and os.path.isfile(self.db_ok):
            return
        
        # Format command
        cmd = [self.mmseqsExe, "createdb", self.fasta, self.db]
        
        # Run DB generation
        print("# MMseqs2 DB generation with: " + " ".join(cmd))
        run_makedb = subprocess.Popen(" ".join(cmd), shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != "":
            raise Exception(f"MMseqs2 database generation of '{self.fasta}' failed with error:\n" +
                            makedberr.decode("utf-8"))
        
        # Create successful completion output flag
        touch_ok(self.db_ok)
    
    def index(self):
        '''
        Indexes as sequence DB for MMSeqs2 if it has not already been done.
        '''
        # Skip if index already exists
        if os.path.isfile(self.index_ok):
            return
        
        # Format command
        cmd = [self.mmseqsExe, "createindex", self.db, self.tmpDir, "--threads", str(self.threads)]
        if self.molecule == "protein":
            cmd += ["--search-type", "1"]
        else:
            cmd += ["--search-type", "3"]
        
        # Run query index
        print("# MMseqs2 DB indexing with: " + " ".join(cmd))
        run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        if indexerr.decode("utf-8") != "":
                raise Exception(f"MMseqs2 indexing of '{self.db}' failed with error:\n" + 
                                indexerr.decode("utf-8"))
        
        # Create successful completion output flag
        touch_ok(self.index_ok)
    
    def clean_all(self):
        '''
        Function to invoke after using a MMSeqs2 database which is no longer
        wanted. It should clean up all files with _seqDB* suffix.
        '''
        for file in os.listdir(self.db_dirname):
            if file.startswith(self.db_basename):
                os.unlink(os.path.join(self.db_dirname, file))

class MMseqs:
    '''
    The MMseqs Class provides access to the MMseqs2 search functionality. It cooperates
    with MM_DB and BLAST_Results to provide a BLAST-like interface for searching a
    database with a query.
    
    Attributes:
        query (REQUIRED) -- a MM_DB object.
        target (REQUIRED) -- a MM_DB object.
        mmseqsDir (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                is found.
        tmpDir (REQUIRED) -- a string location for where MMSeqs2 should keep temp files.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing MMseqs2 search. Defaults to 1.
        # Behavioural parameters
        evalue (OPTIONAL) -- a float or integer value specifying the E-value cut-off to enforce
                             when obtaining MMseqs2 results. Defaults to 1e-5.
        iterations (OPTIONAL) -- an integer >= 1 indicating how many iterations of the search
                                 algorithm should be performed. Defaults to 1.
        sensitivity (OPTIONAL) -- a value in the list [1,2,3,4,5,5.7,6,7,7.5] indicating the
                                  sensitivity of the search. Defaults to 5.7.
        alt_ali (OPTIONAL) -- an integer >= 0 indicating how many alternate alighments should
                              be shown for a query; similar behaviour to BLAST's HSPs. Defaults to 0.
    '''
    def __init__(self, query, target, mmseqsDir, tmpDir):
        self.query = query
        self.target = target
        self.mmseqsDir = mmseqsDir
        self.tmpDir = tmpDir
        
        # Set default attributes
        self.evalue = 1e-5
        self.threads = 1
        self.iterations = 1
        self.sensitivity = 5.7
        self.alt_ali = 0
        
        # Set helper attributes
        self.searchResult = None
        self.isMMseqs = True
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        assert type(value).__name__ == "MM_DB" \
            or hasattr(value, "isMM_DB") and value.isMM_DB == True
        
        self._query = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        assert type(value).__name__ == "MM_DB" \
            or hasattr(value, "isMM_DB") and value.isMM_DB == True
        
        self._target = value
    
    @property
    def mmseqsDir(self):
        return self._mmseqsDir
    
    @mmseqsDir.setter
    def mmseqsDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("mmseqsDir does not point to an existing directory"))
        
        self._mmseqsDir = os.path.abspath(value)
        self.mmseqsExe = os.path.join(self._mmseqsDir, "mmseqs")
    
    @property
    def mmseqsExe(self):
        return self._mmseqsExe
    
    @mmseqsExe.setter
    def mmseqsExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("mmseqsExe does not point to an existing file"))
        
        self._mmseqsExe = value
    
    @property
    def tmpDir(self):
        return self._tmpDir
    
    @tmpDir.setter
    def tmpDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        value = os.path.abspath(value)
        
        if not os.path.isdir(os.path.dirname(value)):
            raise Exception((f"tmpDir's parent location ('{os.path.dirname(value)}') " +
                             "does not exist"))
        
        if not os.path.isdir(value):
            os.mkdir(value)
        
        self._tmpDir = value
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value, "evalue must be >= 0"
        
        self._evalue = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    @property
    def iterations(self):
        return self._iterations
    
    @iterations.setter
    def iterations(self, value):
        assert isinstance(value, int)
        assert 0 < value, "iterations must be a positive integer"
        
        self._iterations = value
    
    @property
    def sensitivity(self):
        return self._sensitivity
    
    @sensitivity.setter
    def sensitivity(self, value):
        assert value in [1,2,3,4,5,5.7,6,7,7.5], \
            "sensitivity must be a value in the list [1,2,3,4,5,5.7,6,7,7.5]"
        
        self._sensitivity = value
    
    @property
    def alt_ali(self):
        return self._alt_ali
    
    @alt_ali.setter
    def alt_ali(self, value):
        assert isinstance(value, int)
        assert 0 <= value, "alt_ali must be >= 0"
        
        self._alt_ali = value
    
    def setup(self):
        '''
        Generates databases and indexes them (if relevant). If the files already
        exist, this function will not take any action.
        '''
        self.query.generate()
        self.target.generate()
        
        self.query.index()
        self.target.index()
    
    def mmseqs(self, outFile, force=False):
        '''
        Performs the MMseqs2 search operation.
        
        Parameters:
            outFile -- a string indicating the location to write the output file to.
            force -- a Boolean indicating whether to force the search to run even if the output file
                     already exists. Defaults to False.
        '''
        self.setup()
        
        # Specify file locations
        outFile = os.path.abspath(outFile) # make sure it's the abspath
        tmpSearchPrefix = outFile + "_search"
        tmpSearchOK = tmpSearchPrefix + ".ok"
        
        # Raise exception if output file already exists
        if os.path.exists(outFile) and os.path.isfile(outFile + ".ok") and force == False:
            raise FileExistsError(f"Output file '{outFile}' already exists!")
        
        # Derive our search type
        if self.query.molecule == "nucleotide" and self.target.molecule == "nucleotide":
            searchType = 3 # mmseqs needs explicit search type for nucleotide>nucleotide searches
        else:
            searchType = 0
        
        # Format search command
        cmd = [
            self.mmseqsExe, "search", self.query.db, self.target.db, tmpSearchPrefix, self.tmpDir,
            "--threads", str(self.threads), "-e", str(self.evalue),
            "--num-iterations", str(self.iterations), "-s", str(self.sensitivity),
            "--alt-ali", str(self.alt_ali), "--search-type", str(searchType)
        ]
        
        # Run MMseqs search
        if not os.path.isfile(tmpSearchOK):
            print("# MMseqs2 searching with: " + " ".join(cmd))
            run_search = subprocess.Popen(" ".join(cmd), shell = True,
                                          stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
            searchout, searcherr = run_search.communicate()
            if searcherr.decode("utf-8") != '':
                raise Exception("MMseqs2 search produced error:\n" + 
                                searcherr.decode("utf-8"))
            touch_ok(tmpSearchOK)
        
        # Format tabulation command
        cmd = [
            self.mmseqsExe, "convertalis", self.query.db, self.target.db, tmpSearchPrefix,
            outFile, "--threads", str(self.threads)
        ]
        
        # Run MMseqs2 convertalis
        print("# Generating MMseqs2 tabular output with: " + " ".join(cmd))
        run_table = subprocess.Popen(" ".join(cmd), shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        tableout, tablerr = run_table.communicate()
        if tablerr.decode("utf-8") != '':
            raise Exception("MMseqs2 convertalis produced error:\n" + 
                            tablerr.decode("utf-8"))
        touch_ok(outFile)
        
        # Clean up temporary search outputs
        outDir = os.path.dirname(outFile)
        for file in os.listdir(outDir):
            if file.startswith(os.path.basename(tmpSearchPrefix)):
                os.unlink(os.path.join(outDir, file))
        
        self.searchResult = outFile
    
    def parse_result(self, searchResultFile=None):
        if searchResultFile == None:
            searchResultFile = self.searchResult
        if searchResultFile == None:
            raise FileNotFoundError(("No search result file specified, and none found from previous " +
                                     "use of this object"))
        
        resultObj = BLAST_Results(searchResultFile)
        resultObj.evalue = self.evalue
        resultObj.parse()
        return resultObj.results

class MM_Clust:
    '''
    The MM_Clust Class behaves as an abstract class for inheritance of basic
    validation logic of parameters shared in common.
    
    Attributes:
        mmDB (REQUIRED) -- a MM_DB object from ZS_BlastIO.
        evalue (OPTIONAL) -- a positive float with a minimum of 0.0 controlling the
                             E-value threshold for clustering; default == 1e-3.
        identity (OPTIONAL) -- a positive float in the range 0.0 -> 1.0 controlling the
                               sequence identity threshold for clustering; default == 0.9.
        cov_pct (OPTIONAL) -- a positive float in the range 0.0 -> 1.0 controlling the
                              amount of aligned residues in both shorter and longer
                              sequences; default == 0.8.
        clust_mode (OPTIONAL) -- a string in the list ["set-cover", "connected-component",
                                 "greedy"], corresponding to modes 0, 1, and 2,3;
                                 default == "set-cover", but optimal is probably
                                 "connected-component".
        threads (OPTIONAL) -- a positive integer for how many threads to use when running
                              MMseqs2 clustering (default==1).
        tmpDir (OPTIONAL) -- a string location for where MMseqs2 should keep temp files;
                             if unspecified, it will use the same location as
                             mmDB.tmpDir.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9, cov_pct=0.8,
                 clust_mode="set-cover", threads=1, tmpDir=None):
        self.mmDB = mmDB
        self.evalue = evalue
        self.identity = identity
        self.cov_pct = cov_pct
        self.clust_mode = clust_mode
        self.threads = threads
        self.tmpDir = tmpDir
    
    @property
    def mmDB(self):
        return self._mmDB
    
    @mmDB.setter
    def mmDB(self, value):
        assert type(value).__name__ == "MM_DB" \
            or hasattr(value, "isMM_DB") and value.isMM_DB == True
        
        self._mmDB = value
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value, "evalue must be >= 0"
        
        self._evalue = value
    
    @property
    def identity(self):
        return self._identity
    
    @identity.setter
    def identity(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value <= 1.0, "identity must be in the range of 0.0 -> 1.0 (inclusive)"
        
        self._identity = value
    
    @property
    def cov_pct(self):
        return self._cov_pct
    
    @cov_pct.setter
    def cov_pct(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value <= 1.0, "cov_pct must be in the range of 0.0 -> 1.0 (inclusive)"
        
        self._cov_pct = value
    
    @property
    def clust_mode(self):
        return self._clust_mode
    
    @clust_mode.setter
    def clust_mode(self, value):
        assert isinstance(value, str)
        assert value in ["set-cover", "connected-component", "greedy"], \
            "clust_mode must be one of the supported Linclust algorithms"
        
        self._clust_mode = ["set-cover", "connected-component", "greedy"].index(value)
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    @property
    def tmpDir(self):
        return self._tmpDir
    
    @tmpDir.setter
    def tmpDir(self, value):
        if value == None:
            value = self.mmDB.tmpDir
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        value = os.path.abspath(value)
        
        if not os.path.isdir(os.path.dirname(value)):
            raise Exception((f"tmpDir's parent location ('{os.path.dirname(value)}') " +
                            "does not exist"))
        
        if not os.path.isdir(value):
            os.mkdir(value)
        
        self._tmpDir = value
    
    def cluster(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def tabulate(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def clean_all(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def parse_tsv(self, tabulatedFile):
        '''
        Parses the output of .tabulate() (aka mmseqs createtsv) and produces
        a cluster dictionary.
        
        Parameters:
            tabulatedFile -- a string indicating the file location of the
                             output from createtsv.
        Returns:
            clusterDict -- a dictionary with structure like:
                           {
                               0: [seqid1, seqid2, ...],
                               1: [ ... ],
                               ...
                           }
        '''
        clusterDict = {}
        clusterNum = -1
        lastCluster = None
        with open(tabulatedFile, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n ")
                if l != "":
                    clustID, seqID = l.split("\t")
                    if clustID != lastCluster:
                        clusterNum += 1
                    lastCluster = clustID
                    
                    clusterDict.setdefault(clusterNum, [])
                    clusterDict[clusterNum].append(seqID)
        return clusterDict

class MM_Linclust(MM_Clust):
    '''
    The MM_Linclust Class provides the logic for running Linclust with an MM_DB object
    as input. The MMseqs executable location and tmpDir will be pulled from the
    mmDB input; a new tmpDir for running linclust can be specified if desired.
    
    Attributes:
        super class attributes -- see docstring for MM_Clust.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9,
                 cov_pct=0.8, clust_mode="set-cover", threads=1, tmpDir=None):
        
        super().__init__(mmDB, evalue, identity, cov_pct,
                        clust_mode, threads, tmpDir)
    
    @property
    def paramsHash(self):
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode)
        return sha256(bytes(strForHash, "utf-8")).hexdigest()
    
    @property
    def db(self):
        return os.path.abspath(f"{self.mmDB.fasta}_{self.paramsHash}_linclustDB")
    
    @property
    def db_basename(self):
        return os.path.basename(self.db)
    
    @property
    def db_dirname(self):
        return os.path.dirname(self.db)
    
    @property
    def db_ok(self):
        return self.db + ".ok"
    
    def cluster(self):
        # Run DB generation & indexing if relevant
        self.mmDB.generate()
        self.mmDB.index()
        
        # Skip if db already exists
        if os.path.isfile(self.db) and os.path.isfile(self.db_ok):
            return
        
        # Format command
        cmd = [
            self.mmDB.mmseqsExe, "linclust", self.mmDB.db, self.db, self.tmpDir,
            "--min-seq-id", str(self.identity), "-c", str(self.cov_pct),
            "-e", str(self.evalue), "--cluster-mode", str(self.clust_mode),
            "--threads", str(self.threads)
        ]
        
        # Clustering
        print("# Running linclust with: " + " ".join(cmd))
        run_linclust = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        linclustout, linclusterr = run_linclust.communicate()
        if linclusterr.decode("utf-8") != '':
            raise Exception('Linclust error text below\n' +
                            linclusterr.decode("utf-8"))
        touch_ok(self.db_ok)
    
    def tabulate(self, outputFileName):
        '''
        Tabulates a linclust database file.
        
        Parameters:
            outputFileName -- a string indicating the file name to write TSV formatted
                              clustering results to.
        '''
        # Skip if table already exists
        if os.path.isfile(outputFileName) and os.path.isfile(outputFileName + ".ok"):
            return
        
        # Format command
        cmd = [self.mmDB.mmseqsExe, "createtsv", self.mmDB.db, self.mmDB.db, self.db, outputFileName]
        
        # Tabulation
        print("# Running table generation with: " + " ".join(cmd))
        run_tabulate = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        tableout, tableerr = run_tabulate.communicate()
        if tableerr.decode("utf-8") != '':
            raise Exception('Linclust tabulation text below\n' +
                            tableerr.decode("utf-8"))
        touch_ok(outputFileName)
    
    def clean_all(self):
        '''
        Function to invoke after performing Linclust, the results of which
        are no longer wanted. It should clean up all files with _linclustDB* suffix.
        '''
        for file in os.listdir(self.db_dirname):
            if file.startswith(self.db_basename):
                os.unlink(os.path.join(self.db_dirname, file))

class MM_Cascade(MM_Clust):
    '''
    The MM_Cascade Class provides the logic for running cascaded MMseqs2 clustering
    with an MM_DB object as input. The MMseqs executable location and tmpDir will be pulled
    from the mmDB input; a new tmpDir for running cascaded clusering can be specified if desired.
    
    Attributes:
        super class attributes -- see docstring for MM_Clust.
        sensitivity -- a float in the list [1,2,3,4,5,5.7,6,7,7.5]; default == 4.0,
                       but recommended to use 5.7 or 7.5.
        cluster_steps -- an int minimally bounded at 1 for how many cascaded clustering
                         steps to run; default == 3, recommended to keep that.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9,
                 cov_pct=0.8, clust_mode="set-cover",
                 threads=1, tmpDir=None,
                 sensitivity=4.0, cluster_steps=3):
        
        super().__init__(mmDB, evalue, identity, cov_pct,
                        clust_mode, threads, tmpDir)
        
        self.sensitivity = sensitivity
        self.cluster_steps = cluster_steps
    
    @property
    def sensitivity(self):
        return self._sensitivity
    
    @sensitivity.setter
    def sensitivity(self, value):
        assert value in [1,2,3,4,5,5.7,6,7,7.5], \
            "sensitivity must be a value in the list [1,2,3,4,5,5.7,6,7,7.5]"
        
        self._sensitivity = value
    
    @property
    def cluster_steps(self):
        return self._cluster_steps
    
    @cluster_steps.setter
    def cluster_steps(self, value):
        assert isinstance(value, int)
        assert 0 < value, "cluster_steps must be a positive integer"
        
        self._cluster_steps = value
    
    @property
    def paramsHash(self):
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode) + str(self.sensitivity)
        return sha256(bytes(strForHash, "utf-8")).hexdigest()
    
    @property
    def db(self):
        return os.path.abspath(f"{self.mmDB.fasta}_{self.paramsHash}_clustDB")
    
    @property
    def db_basename(self):
        return os.path.basename(self.db)
    
    @property
    def db_dirname(self):
        return os.path.dirname(self.db)
    
    @property
    def db_ok(self):
        return self.db + ".ok"
    
    def cluster(self):
        # Run DB generation & indexing if relevant
        self.mmDB.generate()
        self.mmDB.index()
        
        # Skip if db already exists
        if os.path.isfile(self.db) and os.path.isfile(self.db_ok):
            return
        
        # Format command
        cmd = [
            self.mmDB.mmseqsExe, "cluster", self.mmDB.db, self.db, self.tmpDir,
            "--min-seq-id", str(self.identity), "-c", str(self.cov_pct),
            "-e", str(self.evalue), "--cluster-mode", str(self.clust_mode),
            "-s", str(self.sensitivity), "--cluster-steps", str(self.cluster_steps),
            "--threads", str(self.threads)
        ]
        
        # Clustering
        print("# Running cascaded clustering with: " + " ".join(cmd))
        run_cluster = subprocess.Popen(" ".join(cmd), shell = True,
                                       stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        clustout, clusterr = run_cluster.communicate()
        if clusterr.decode("utf-8") != '':
            raise Exception('Cascaded clustering error text below\n' +
                            clusterr.decode("utf-8"))
        touch_ok(self.db_ok)
    
    def tabulate(self, outputFileName):
        '''
        Tabulates a cascaded database file.
        
        Parameters:
            outputFileName -- a string indicating the file name to write TSV formatted
                              clustering results to.
        '''
        # Skip if table already exists
        if os.path.isfile(outputFileName) and os.path.isfile(outputFileName + ".ok"):
            return
        
        # Format command
        cmd = [self.mmDB.mmseqsExe, "createtsv", self.mmDB.db, self.mmDB.db, self.db, outputFileName]
        
        # Tabulation
        print("# Running table generation with: " + " ".join(cmd))
        run_tabulate = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        tableout, tableerr = run_tabulate.communicate()
        if tableerr.decode("utf-8") != '':
            raise Exception('Cascaded clustering tabulation text below\n' +
                            tableerr.decode("utf-8"))
        touch_ok(outputFileName)
    
    def clean_all(self):
        '''
        Function to invoke after performing cascaded clustering, the results of which
        are no longer wanted. It should clean up all files with _clustDB* suffix.
        '''
        for file in os.listdir(self.db_dirname):
            if file.startswith(self.db_basename):
                os.unlink(os.path.join(self.db_dirname, file))
