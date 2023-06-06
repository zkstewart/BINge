import os

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
