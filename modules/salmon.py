import os, subprocess, platform

class Salmon:
    '''
    The Salmon Class provides access to the salmon read mapping functionality. It cooperates
    with a Salmon_DB class object to perform the mapping.
    
    Attributes:
        reads (REQUIRED) -- a list containing one string (single-end) or two strings (paired-end)
                            indicating the location of FASTQ file(s) to map.
        target (REQUIRED) -- a Salmon_DB object.
        salmonDir (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                is found.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing MMseqs2 search. Defaults to 1.
    '''
    def __init__(self, reads, target, salmonDir, threads=1):
        self.reads = reads
        self.target = target
        self.salmonDir = salmonDir
        self.threads = threads
        
        # Set helper attributes
        self.isSetup = False
        self.searchResult = None
        self.isSalmon = True
    
    @property
    def reads(self):
        return self._reads
    
    @reads.setter
    def reads(self, value):
        assert isinstance(value, list)
        assert len(value) == 1 or len(value) == 2, \
            "reads must be a list of one or two strings"
        for read in value:
            assert isinstance(read, str)
            if not os.path.isfile(read):
                raise FileNotFoundError(f"Read file does not exist: {read}")
        
        self._reads = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        assert type(value).__name__ == "Salmon_DB" \
            or type(value).__name__ == "ZS_MapIO.Salmon_DB" \
            or hasattr(value, "isSalmon_DB") and value.isSalmon_DB == True
        
        self._target = value
    
    @property
    def salmonDir(self):
        return self._salmonDir
    
    @salmonDir.setter
    def salmonDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("salmonDir does not point to an existing directory"))
        
        self._salmonDir = value
        self.salmonExe = os.path.join(value, "salmon")
    
    @property
    def salmonExe(self):
        return self._salmonExe
    
    @salmonExe.setter
    def salmonExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("salmonExe does not point to an existing file"))
        
        self._salmonExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    def setup(self):
        '''
        Indexes database (if relevant)
        '''
        if not self.target.index_exists():
            self.target.index()
        
        self.isSetup = True
    
    def quant(self, outputDirectory):
        '''
        Performs the salmon quant operation.
        
        Parameters:
            outputDirectory -- a string indicating the location to
                               write outputs to.
        '''
        if not self.isSetup:
            self.setup()
        
        # Specify file locations
        reads = [os.path.abspath(read) for read in self.reads]
        target = os.path.abspath(f"{self.target.fasta}.salmonDB")
        
        # Format command
        cmd = [
            self.salmonExe, "quant", "--threads", str(self.threads),
            "--index", target, "--libType", "A", "-o", outputDirectory
        ]
        
        # Add read files
        if len(reads) == 1:
            cmd += ["-r", reads[0]]
        elif len(reads) == 2:
            cmd += ["-1", reads[0], "-2", reads[1]]
        else:
            raise Exception("Incorrect number of read files provided")
        
        # Perform salmon search
        if platform.system() != "Windows":
            run_salmon = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_salmon = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        
        salmonout, salmonerr = run_salmon.communicate()
        if "writing output" not in salmonerr.decode("utf-8"):
            raise Exception('salmon searching error text below\n' +
                            salmonerr.decode("utf-8"))

class Salmon_DB:
    '''
    The Salmon_DB Class encapsulates the logic of indexing a FASTA file using salmon.
    
    Attributes:
        fasta (REQUIRED) -- a string indicating the location of a FASTA file.
        salmonDir (REQUIRED) -- a string indicating the location of salmon binaries.
    '''
    def __init__(self, fasta, salmonDir, threads=1):
        self.fasta = fasta
        self.salmonDir = salmonDir
        self.threads = threads
        
        self.isSalmon_DB = True # flag to check object type
    
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
    def salmonDir(self):
        return self._salmonDir
    
    @salmonDir.setter
    def salmonDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("salmonDir does not point to an existing directory"))
        
        self._salmonDir = value
        self.salmonExe = os.path.join(value, "salmon")
    
    @property
    def salmonExe(self):
        return self._salmonExe
    
    @salmonExe.setter
    def salmonExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("salmonExe does not point to an existing file"))
        
        self._salmonExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("Number of threads must be more than 0")
        
        self._threads = value
    
    def index_exists(self):
        '''
        Checks for the expected index directory.
        
        Returns:
            dbExists -- a Boolean where True means the database exists,
                        and False means it does not exist.
        '''
        return os.path.isdir(f"{self.fasta}.salmonDB")
    
    def index(self):
        '''
        Makes a database out of the .fasta value for use in salmon mapping.
        '''
        # Skip if index exists
        if self.index_exists():
            raise FileExistsError("salmon index already exists!")
        
        # Specify file locations
        fasta = os.path.abspath(self.fasta)
        
        # Format command
        cmd = [
            self.salmonExe, "index", "--threads", str(self.threads),
            "--transcripts", fasta, "--index", f"{fasta}.salmonDB"
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
        if "done building index" not in indexerr.decode("utf-8"):
            raise Exception('salmon indexing error text below\n' + indexerr.decode("utf-8"))

class EquivalenceClassCollection():
    def __init__(self):
        self.ids = {} # transcript IDs to EC index dict
        self.samples = [] # sample names
        self.ec = {} # equivalence class dicts
        
        self.numTranscripts = None # for checking EC file compatibility
    
    def parse_eq_file(self, eqFile, sample):
        # Validate input parameters
        assert os.path.isfile(eqFile), \
            f"Cannot parse '{eqFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Figure out if we should hold onto transcript IDs or not
        if self.ids == {}:
            needsIDs = True
        else:
            needsIDs = False
        
        # Parse the file proper
        with open(eqFile, "r") as fileIn:
            
            # Parse the first two lines of the file
            numTranscripts = int(fileIn.readline().rstrip("\r\n "))
            numECs = int(fileIn.readline().rstrip("\r\n "))
            
            # Check compatibility of EC files
            if self.numTranscripts == None:
                self.numTranscripts = numTranscripts
            else:
                assert self.numTranscripts == numTranscripts, \
                    "EC files are incompatible since transcript numbers differ!"
            
            # Loop through transcript ID lines
            for i in range(numTranscripts):
                transcriptID = fileIn.readline().rstrip("\r\n ")
                if needsIDs:
                    self.ids[transcriptID] = i
            
            # Loop through EC lines
            thisEC = { i:0.0 for i in range(0, numTranscripts) } # holds onto counts per transcript for this file
            for _ in range(numECs):
                sl = list(map(int, fileIn.readline().rstrip("\r\n ").split("\t")))
                numIDs, idIndices, numReads = sl[0], sl[1:-1], sl[-1]
                numReads = numReads / numIDs
                
                # Store relevant data
                for idIndex in idIndices:
                    thisEC[idIndex] += numReads
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store EC counts inside parent structure
        for ecIndex, ecCount in thisEC.items():
            self.ec.setdefault(ecIndex, [])
            self.ec[ecIndex].append(ecCount)
    
    def get_transcript_count(self, ecValue):
        '''
        Parameters:
            ecValue -- an integer identifying the equivalence class index
                       you want to get counts for, OR a string identifying
                       the transcript ID you want to get counts for.
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        if isinstance(ecValue, int):
            return self.ec[ecValue]
        else:
            return self.ec[self.ids[ecValue]]
    
    def __repr__(self):
        return "<EquivalenceClassCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

class QuantCollection():
    def __init__(self):
        self.samples = [] # sample names
        self.quant = {} # quant count dicts
        
        self.numTranscripts = None # for checking quant file compatibility
    
    def parse_quant_file(self, quantFile, sample):
        # Validate input parameters
        assert os.path.isfile(quantFile), \
            f"Cannot parse '{quantFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Parse the file proper
        firstLine = True
        thisQuant = {}
        with open(quantFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                
                # Handle header line
                if firstLine == True:
                    assert sl == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"], \
                        "Quant file appears to lack the expected header? Cannot parse."
                    firstLine = False
                # Handle content lines
                else:
                    name, _, _, _, numReads = sl # might error here if file format is bad
                    thisQuant[name] = float(numReads)
        
        # Check compatibility of quant file with existing ones
        if self.numTranscripts == None:
            self.numTranscripts = len(thisQuant)
        else:
            assert self.numTranscripts == len(thisQuant), \
                "Quant files are incompatible since transcript numbers differ!"
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store counts inside parent structure
        for name, numReads in thisQuant.items():
            self.quant.setdefault(name, [])
            self.quant[name].append(float(numReads))
    
    def get_transcript_count(self, seqID):
        '''
        Parameters:
            seqID -- a string identifying the transcript ID you want to get counts for.
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        return self.quant[seqID]
    
    def __repr__(self):
        return "<QuantCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

class DGEQuantCollection():
    '''
    Similar to QuantCollection, but it stores more data which is necessary for DGE
    analysis i.e., the effective length and TPM values.
    '''
    def __init__(self):
        self.samples = [] # sample names
        self.quant = {}
        self.length = {}
        self.tpm = {}
        
        self.numTranscripts = None # for checking quant file compatibility
    
    def parse_quant_file(self, quantFile, sample):
        # Validate input parameters
        assert os.path.isfile(quantFile), \
            f"Cannot parse '{quantFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Parse the file proper
        firstLine = True
        thisQuant = {}
        
        with open(quantFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                
                # Handle header line
                if firstLine == True:
                    assert sl == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"], \
                        "Quant file appears to lack the expected header? Cannot parse."
                    firstLine = False
                # Handle content lines
                else:
                    name, _, effectiveLength, tpm, numReads = sl # might error here if file format is bad
                    thisQuant[name] = [float(numReads), float(effectiveLength), float(tpm)]
        
        # Check compatibility of quant file with existing ones
        if self.numTranscripts == None:
            self.numTranscripts = len(thisQuant)
        else:
            assert self.numTranscripts == len(thisQuant), \
                "Quant files are incompatible since transcript numbers differ!"
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)
        
        # Store counts inside parent structure
        for name, values in thisQuant.items():
            self.quant.setdefault(name, [])
            self.length.setdefault(name, [])
            self.tpm.setdefault(name, [])
            
            numReads, effectiveLength, tpm = values
            self.quant[name].append(numReads)
            self.length[name].append(effectiveLength)
            self.tpm[name].append(tpm)
    
    def get_transcript_count(self, seqID):
        '''
        Parameters:
            seqID -- a string identifying the transcript ID you want to get counts for.
        Returns:
            transcriptCount -- a list containing each sample's count for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        return self.quant[seqID]
    
    def get_transcript_effective_length(self, seqID):
        '''
        Parameters:
            seqID -- a string identifying the transcript ID you want to get counts for.
        Returns:
            effectiveLength -- a list containing each sample's effective length for
                               the given transcript; order is equivalent to
                               self.samples.
        '''
        return self.length[seqID]
    
    def get_transcript_tpm(self, seqID):
        '''
        Parameters:
            seqID -- a string identifying the transcript ID you want to get counts for.
        Returns:
            get_transcript_tpm -- a list containing each sample's TPM abundance for
                                  the given transcript; order is equivalent to
                                  self.samples.
        '''
        return self.tpm[seqID]
    
    def __repr__(self):
        return "<DGEQuantCollection object;num_samples={0};num_transcripts={1}>".format(
            len(self.samples), self.numTranscripts
        )

class SalmonQC():
    '''
    Parses salmon output directories for log files to extract QC information.
    '''
    def __init__(self, salmonDirs):
        self.numReads = {} # read counts from quant files
        self.totalReads = {} # total read counts from quant files
        self.mapPct = {} # mapping percentages from log files
        
        self.salmonDirs = salmonDirs
    
    @property
    def salmonDirs(self):
        return self._salmonDirs
    
    @salmonDirs.setter
    def salmonDirs(self, salmonDirs):
        # Validate input parameters
        assert isinstance(salmonDirs, list), \
            f"salmonDirs must be a list of directory paths!"
        
        # Parse the files
        for salmonDir in salmonDirs:
            self.parse_salmon_dir(salmonDir)
    
    def parse_quant_file(self, quantFile, sample):
        # Validate input parameters
        assert os.path.isfile(quantFile), \
            f"Cannot parse '{quantFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.numReads, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        self.numReads[sample] = {}
        
        # Parse the file
        firstLine = True
        total = 0
        with open(quantFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                
                if firstLine == True:
                    firstLine = False
                else:
                    name, _, _, _, numReads = sl
                    numReads = float(numReads)
                    self.numReads[sample][name] = numReads
                    total += numReads
        self.totalReads[sample] = total
    
    def parse_log_file(self, logFile, sample):
        # Validate input parameters
        assert os.path.isfile(logFile), \
            f"Cannot parse '{logFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.mapPct, \
            f"Sample name must be a string value that uniquely identifies this sample!"
        
        # Parse the file
        mappingRate = None
        with open(logFile, "r") as fileIn:
            for line in fileIn:
                if "[info] Mapping rate =" in line:
                    mappingRate = line.rstrip("\r\n ").split(" = ")[1].rstrip("%")
                    foundMappingRate = True
                    break
        if mappingRate == None:
            raise ValueError(f"Cannot locate mapping rate in '{logFile}'!")
        
        # Store the mapping rate
        self.mapPct[sample] = float(mappingRate)
    
    def parse_salmon_dir(self, salmonDir):
        # Validate input parameters
        assert os.path.isdir(salmonDir), \
            f"Cannot parse '{salmonDir}' as directory does not exist!"
        
        # Derive sample name from directory
        sample = os.path.basename(salmonDir)
        
        # Locate quant and log files
        quantFile = os.path.join(salmonDir, "quant.sf")
        if not os.path.isfile(quantFile):
            raise FileNotFoundError(f"Cannot locate quant file in '{salmonDir}'!")
        
        logFile = os.path.join(salmonDir, "logs", "salmon_quant.log")
        if not os.path.isfile(logFile):
            raise FileNotFoundError(f"Cannot locate log file in '{salmonDir}'!")
        
        # Parse quant and log files
        self.parse_quant_file(quantFile, sample)
        self.parse_log_file(logFile, sample)
    
    def get_clustered_qc(self, clusterDict):
        '''
        Computes the percentage of reads that are represented in the clustered sequences
        for each sample, and then adjusts the total mapping percentage on this basis.
        
        Parameters:
            clusterDict -- a dictionary with format like:
                           { 0: [ [seqid1, seqid2, ...], "binned" ],
                             3: [ [...], "unbinned"],
                             12: ...,
                             ...
                           }
        '''
        adjustedQC = {}
        for sample, numReads in self.totalReads.items():
            # Get the number of reads that clustered for this sample
            clusteredNumReads = 0
            for clusterNum, seqIDs in clusterDict.items():
                for seqID in seqIDs:
                    if seqID in self.numReads[sample]:
                        clusteredNumReads += self.numReads[sample][seqID]
                    else:
                        continue # salmon can sometimes omit transcripts from quant.sf files
            
            # Get the percentage of reads that clustered out of the amount that mapped initially
            retainedPct = (clusteredNumReads / numReads) * 100
            
            # Adjust the mapped percentage according to the original mapping percentage
            clusteredPct = self.mapPct[sample] * (retainedPct / 100)
            
            # Store the adjusted QC value
            adjustedQC[sample] = [retainedPct, clusteredPct]
        
        return adjustedQC
    
    def __repr__(self):
        return "<SalmonQC object;num_samples={0}>".format(
            len(self.mapPct)
        )
