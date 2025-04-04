import os, shutil, subprocess, platform

class CDHIT:
    '''
    The CDHIT Class provides easy access to CD-HIT functionality to perform
    redundancy reduction of FASTA files. These can be in the form of providing
    the location of a FASTA file as string, or with ZS_SeqIO modules i.e.,
    FASTA and FastASeq.
    
    Attributes:
        fasta (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        molecule (REQUIRED) -- a string, indicating whether the FASTA file contains protein
                               or nucleotide sequences. Value must be in the list:
                               [n, nucl, nucleotide, p, prot, protein]
        cdhitDir (OPTIONAL) -- a string indicating the location of the CD-HIT executables,
                               if they're not available from the PATH variable.
    '''
    def __init__(self, fasta, molecule, cdhitDir=None):
        # Validate input types and values
        assert type(fasta).__name__ == "str"
        
        assert isinstance(molecule, str)
        assert molecule.lower() in ["n", "nucl", "nucleotide", "p", "prot", "protein"]
        
        assert isinstance(cdhitDir, str) or cdhitDir == None
        
        # Validate that fasta exists if specified as a string (file location)
        if type(fasta).__name__ == "str" and not os.path.isfile(fasta):
            raise Exception("fasta parameter is a string, but does not point to an existing file location")
        self.fasta = fasta
        
        # Validate that cdhitDir exists and contains the relevant executables if specified as a string
        if isinstance(cdhitDir, str):
            assert os.path.isdir(cdhitDir), "cdhitDir does not exist"
            assert os.path.isfile(os.path.join(cdhitDir, "cd-hit")) or os.path.isfile(os.path.join(cdhitDir, "cd-hit.exe")), "cd-hit executable not found at '{0}'".format(cdhitDir)
            assert os.path.isfile(os.path.join(cdhitDir, "cd-hit-est")) or os.path.isfile(os.path.join(cdhitDir, "cd-hit-est.exe")), "cd-hit-est executable not found at '{0}'".format(cdhitDir)
        else:
            assert shutil.which("cd-hit") != None, "cd-hit executable not found in path"
            assert shutil.which("cd-hit-est") != None, "cd-hit-est executable not found in path"
        self.cdhitDir = cdhitDir
        
        # Coerce molecule into a consistent value
        self.molecule = "nucleotide" if molecule.lower() in ["n", "nucl", "nucleotide"] else "protein"
        
        # Set default attributes
        self.identity = 0.9 # implicitly sets self.word_length
        self.local = False
        self.shorter_cov_pct = 0.0
        self.longer_cov_pct = 0.0
        self.mem = 1000 # 1GB by default
        self.threads = 1
        self.clean = True
        self.description_length = 0
        
        # Results storage values
        self.resultFasta = None
        self.resultClusters = None
    
    @property
    def identity(self):
        return self._identity
    
    @identity.setter
    def identity(self, num):
        '''
        Relates to the -c parameter for CD-HIT. Controls the sequence identity
        threshold for clustering sequences together. This method also controls
        the word length setting, since the word length should always correspond
        to the identity threshold anyway.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not 0<=num<=1:
            raise Exception("Identity (-c) must be between 0 and 1 (inclusive)")
        
        self._identity = num
        
        # Set word length intelligently, depending on molecule type
        if self.molecule == "nucleotide":
            if 0.9<=num<=1.0:
                self.word_length = 8
            elif 0.88<=num<=0.9:
                self.word_length = 7
            elif 0.85<=num<=0.88:
                self.word_length = 6
            elif 0.80<=num<=0.85:
                self.word_length = 5
            elif 0.75<=num<=0.80:
                self.word_length = 4
            else:
                self.word_length = 3 # is this okay? I don't know
        else:
            if 0.7<=num<=1.0:
                self.word_length = 5
            elif 0.6<=num<=0.7:
                self.word_length = 4
            elif 0.5<=num<=0.6:
                self.word_length = 3
            else:
                self.word_length = 2 # assuming this is also okay
    
    @property
    def word_length(self):
        return self._word_length
    
    @word_length.setter
    def word_length(self, num):
        '''
        Relates to the -n parameter for CD-HIT. Controls algorithmic behaviour
        of CD-HIT, and should be set to a value according to the identity parameter.
        
        Usually, you should allow the identity setter to automatically specify the
        word length, but this setter allows you to override that behaviour.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 1 and
                   maximally bounded at 11.
        '''
        assert isinstance(num, int)
        if not (1 <= num <= 11):
            raise Exception("word length (-n) must be in the range of 1 to 11 (inclusive)")
        
        self._word_length = num
    
    @property
    def mem(self):
        return self._mem
    
    @mem.setter
    def mem(self, num):
        '''
        Relates to the -M parameter for CD-HIT. Controls how much memory is allowed
        to be used by the algorithm. Number is in megabytes.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 100.
        '''
        assert isinstance(num, int)
        if not num >= 100:
            raise Exception("mem (-M) must be at least 100 megabytes")
        
        self._mem = num
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, num):
        '''
        Relates to the -T parameter for CD-HIT. Controls how many threads CD-HIT
        will use.
        
        0 means it will use as many threads as there are CPU cores.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 0.
        '''
        assert isinstance(num, int)
        if not num >= 0:
            raise Exception("threads (-T) must be greater than or equal to 0")
        
        self._threads = num
    
    @property
    def clean(self):
        return self._clean
    
    @clean.setter
    def clean(self, value):
        '''
        This method allows the clean attribute to be set, which controls
        whether CD-HIT output files are kept or not.
        
        Parameters:
            clean -- a Boolean of True or False
        '''
        assert isinstance(value, bool)
        
        self._clean = value
    
    @property
    def description_length(self):
        return self._description_length
    
    @description_length.setter
    def description_length(self, num):
        '''
        Relates to the -d parameter for CD-HIT. Controls the presentation of
        the output cluster file.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 0.
        '''
        assert isinstance(num, int)
        if not num >= 0:
            raise Exception("description length (-d) must be >= 0")
        
        self._description_length = num
    
    def set_local(self):
        '''
        Relates to the -G parameter for CD-HIT. Changes how identity is scored
        to only consider the local region in which an alignment was found.
        
        This is NOT the default CD-HIT behaviour.
        '''
        self.local = True
    
    def set_global(self):
        '''
        Relates to the -G parameter for CD-HIT. Changes how identity is scored
        to consider the entire sequence length as part of the calculation.
        
        This is the default CD-HIT behaviour
        '''
        self.local = False
    
    def set_shorter_cov_pct(self, num):
        '''
        Relates to the -aS parameter for CD-HIT. Controls minimum threshold
        of sequence that must align well from the shorter sequence. For example,
        a value of 0.5 means at least half of the shorter sequence in a pairwise
        alignment must align against the longer sequence. If it's less, it
        will be discarded.
        
        0.0 means no threshold will be enforced.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not (0<=num<=1):
            raise Exception("shorter_cov_pct (-aS) must be between 0 and 1 (inclusive)")
        
        self.shorter_cov_pct = num
    
    def set_longer_cov_pct(self, num):
        '''
        Relates to the -aL parameter for CD-HIT. Controls minimum threshold
        of sequence that must align well from the longer sequence. For example,
        a value of 0.5 means at least half of the longest sequence in a pairwise
        alignment must align against the shorter sequence. If it's less, it
        will be discarded.
        
        0.0 means no threshold will be enforced.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not (0<=num<=1):
            raise Exception("longer_cov_pct (-aL) must be between 0 and 1 (inclusive)")
        
        self.longer_cov_pct = num

    def cdhit(self, fasta, outputDir, outputFasta):
        '''
        Performs the CD-HIT operation.
        
        This method does not use .fasta because it may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            fasta -- a string indicating the location of a FASTA file to cluster.
            outputDir -- a string indicating the location to write the output FASTA and .clstr files to.
            outFile -- a string indicating the name for the output FASTA file. File must not already exist!
        Returns:
            cmd -- a string indicating the command executed for CD-HIT clustering.
        '''
        # Validate parameters
        assert isinstance(fasta, str)
        assert isinstance(outputDir, str)
        assert isinstance(outputFasta, str)
        
        assert os.path.isfile(fasta), "fasta file does not exist"
        assert os.path.isdir(outputDir), "output directory does not exist"
        assert os.path.basename(outputFasta) == outputFasta, \
            "output fasta file needs to be just the file name; its location is specified in the outputDir method parameter"
        assert not os.path.isfile(os.path.join(outputDir, outputFasta)), \
            f"'{os.path.join(outputDir, outputFasta)}' already exists; cdhit method won't overwrite it"
        
        # Figure out which CD-HIT executable we're using
        if self.molecule == "nucleotide":
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit-est')
        else:
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit')
        
        # Begin formatting cmd, converting to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(fasta)
            outputFile = convert_windows_to_wsl_path(os.path.join(outputDir, outputFasta))
        else:
            outputFile = os.path.join(outputDir, outputFasta)
        cmd = base_subprocess_cmd(program)
        
        # Format cmd and run it
        cmd += list(map(str, ["-i", fasta, "-o", outputFile,
               "-c", self.identity, "-n", self.word_length, "-G", "0" if self.local else "1", 
               "-aS", self.shorter_cov_pct, "-aL", self.longer_cov_pct, "-M", self.mem,
               "-T", self.threads, "-d", self.description_length]))
        
        if platform.system() != "Windows":
            run_cdhit = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_cdhit = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        
        cdout, cderr = run_cdhit.communicate()
        if cderr.decode("utf-8") != '':
            raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))
        
        return cmd
    
    @staticmethod
    def parse_clstr_file(clstrFile):
        '''
        For this to be effective, you should make sure CD-HIT was run with -d 0
        so as to give the sequence ID in a format expected here.
        
        Parameters:
            clstrFile -- a string pointing to the location of a CD-HIT output cluster
                         file.
        Returns:
            clstrDict -- a dictionary with structure like:
                         {
                             0: [seqid1, seqid2, ...],
                             1: [ ... ],
                             ...
                         }
        '''
        clstrDict = {}
        with open(clstrFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split()
                
                # Handle cluster ID lines
                if line.startswith(">"):
                    thisCluster = int(sl[1])
                    clstrDict[thisCluster] = []
                
                # Handle content lines
                else:
                    seqID = sl[2].strip(">.")
                    clstrDict[thisCluster].append(seqID)
        return clstrDict
    
    def get_cdhit_results(self, workingDir=".", returnFASTA=True, returnClusters=False):
        '''
        This function pipelines the process of obtaining CD-HIT results. Intermediate files are
        deleted automatically, and hence this function will only result in the return of the
        clustered FASTA object.
        
        Parameters:
            workingDir -- a string indicating the location to write CD-HIT results to.
            returnFASTA -- a boolean indicating whether we should parse the output FASTA
                           file at the end of this, storing its value in .resultFasta as
                           a ZS_SeqIO.FASTA object.
            returnClusters -- a boolean indicating whether we should parse the .clstr
                              file at the end of this, storing its value in .resultClusters
                              as a dictionary.
        Returns:
            FASTA_obj -- a ZS_SeqIO.FASTA object of the clustered CD-HIT results.
            cdhitResultFile -- a string indicating the file name of the results file. If self.clean is True,
                               this will instead return None.
        '''
        assert os.path.isdir(workingDir), \
            "workingDir must already exist, or just leave it as default to write to current working directory"
        
        # Validate that running this function will result in some sort of output
        if self.clean == True:
            if returnFASTA == False and returnClusters == False:
                raise ValueError(("get_cdhit_results will not return you anything " +
                                  "at the end of running, since .clean is set to True " +
                                  "and you set both returnFASTA and returnClusters to False. " +
                                  "Since there's no point running, I'm not going to."))
        
        # Get file name after data type coercion
        f, fIsTemporary = Conversion.get_filename_for_input_sequences(self.fasta)
        
        # Get hash for temporary file creation
        tmpHash = Conversion.get_hash_for_input_sequences(f)
        
        # Run CD-HIT
        tmpResultName = tmp_file_name_gen("cdhit_result_tmp" + tmpHash, "fasta")
        self.cdhit(f, workingDir, tmpResultName) # "." for working directory being the current one
        
        # Parse CD-HIT results if desired
        if returnFASTA == True:
            result_FASTA_obj = FASTA(tmpResultName)
        else:
            result_FASTA_obj = None
        
        # Clean up f temporary file
        if fIsTemporary:
            os.unlink(f)
        
        # Parse clstr file if desired
        if returnClusters == True:
            clstrDict = CDHIT.parse_clstr_file(tmpResultName + ".clstr")
        else:
            clstrDict = None
        
        # Store results
        self.resultFasta = result_FASTA_obj
        self.resultClusters = clstrDict
        
        # Clean up results (if relevant) and return the output file name
        if self.clean:
            os.unlink(tmpResultName)
            os.unlink(tmpResultName + ".clstr")
            return None
        # Or just return results
        else:
            return tmpResultName
    
    def __repr__(self):
        return (
            f"<CDHIT object;identity={self.identity};local={self.local};" +
            f"shorter_cov_pct={self.shorter_cov_pct};longer_cov_pct={self.longer_cov_pct};" +
            f"mem={self.mem};threads={self.threads};clean={self.clean};" +
            f"description_length={self.description_length};" +
            f"resultFASTA contains data={'NO' if self.resultFasta == None else 'YES'};" +
            f"resultClusters contains data={'NO' if self.resultClusters == None else 'YES'}"
        )
