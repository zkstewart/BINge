import os, re, pickle
from Bio import SeqIO
from pyfaidx import Fasta
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .thread_workers import BasicProcess
from .validation import handle_symlink_change, touch_ok

# Define classes
class FastaParser:
    def __init__(self, file):
        self.file = file
    
    def __iter__(self):
        with open(self.file, 'r') as fastaFile:
            for title, seq in SimpleFastaParser(fastaFile):
                yield title, seq

class FastaCollection:
    '''
    Wrapper for pyfaidx Fasta objects which allows multiple to be combined
    and queried as one logical entity.
    
    Parameters:
        fastaFiles -- a list of strings pointing to the locations of FASTA files
                      which are to be loaded in using pyfaidx.Fasta
    '''
    def __init__(self, fastaFiles):
        self.fastaFiles = fastaFiles
        self.records = []
        self.pipePrefixes = []
        
        self._parse_fastas()
    
    def _parse_fastas(self):
        for fastaFile in self.fastaFiles:
            self.records.append(Fasta(fastaFile))
            # Extract pipe prefixes from FASTA titles
            pipePrefixes = set()
            with open(fastaFile, "r") as fileIn:
                for line in fileIn:
                    if line.startswith(">"):
                        seqPrefix = line[1:].split(None, 1)[0]
                        if "|" in seqPrefix:
                            pipePrefix = seqPrefix.split("|")[0] + "|"
                            pipePrefixes.add(pipePrefix)
            self.pipePrefixes.append(pipePrefixes)
    
    def __getitem__(self, key):
        for i, records in enumerate(self.records):
            try:
                return records[key]
            except:
                for pipePrefix in self.pipePrefixes[i]:
                    try:
                        return records[pipePrefix + key]
                    except:
                        pass
        raise KeyError(f"'{key}' not found in collection")
    
    def __contains__(self, key):
        try:
            self[key] # __getitem__ raises exception if the key isn't found
            return True
        except:
            return False
    
    def __iter__(self):
        for records in self.records:
            yield from records
    
    def __repr__(self):
        return (f"<FastaCollection object;num_records='{len(self.records)}';" +
                f"fastaFiles={self.fastaFiles}"
        )

class ORF_Find:
    '''
    Copy-pasted reimplementation of the biopython_orf_find.py script's logic with
    some minor modifications to remove unneeded functionality.
    '''
    def __init__(self, fastaFile):
        '''
        Params:
            fastaFile -- a string indicating the location of a FASTA file.
        '''
        self.fastaFile = fastaFile
        self.startCodon = re.compile(r'^.*?(M.*)')
        self.xRegex = re.compile(r'X+')
        
        # Set default values
        self._minProLen = 30
        self._maxProLen = 0
        self._hitsToPull = 3
        self._altCodonStringency = 49
        self._noCodonStringency = 99
        self._sequenceType = "prot"
        self._unresolvedCodon = 0
        self._translationTable = 1
    
    @property
    def minProLen(self):
        return self._minProLen
    
    @minProLen.setter
    def minProLen(self, value):
        assert isinstance(value, int)
        assert value > 0, "minProLen must be greater than 0"
        self._minProLen = value
    
    @property
    def maxProLen(self):
        return self._maxProLen
    
    @maxProLen.setter
    def maxProLen(self, value):
        assert isinstance(value, int)
        assert value >= 0, "maxProLen must be greater than or equal to 0"
        self._maxProLen = value
    
    @property
    def hitsToPull(self):
        return self._hitsToPull
    
    @hitsToPull.setter
    def hitsToPull(self, value):
        assert isinstance(value, int)
        assert value > 0, "hitsToPull must be greater than 0"
        self._hitsToPull = value
    
    @property
    def altCodonStringency(self):
        return self._altCodonStringency
    
    @altCodonStringency.setter
    def altCodonStringency(self, value):
        assert isinstance(value, int)
        assert value >= 0, "altCodonStringency must be greater than or equal to 0"
        self._altCodonStringency = value
    
    @property
    def noCodonStringency(self):
        return self._noCodonStringency
    
    @noCodonStringency.setter
    def noCodonStringency(self, value):
        assert isinstance(value, int)
        assert value >= 0, "noCodonStringency must be greater than or equal to 0"
        self._noCodonStringency = value
    
    @property
    def unresolvedCodon(self):
        return self._unresolvedCodon
    
    @unresolvedCodon.setter
    def unresolvedCodon(self, value):
        assert isinstance(value, int)
        assert value >= 0, "unresolvedCodon must be greater than or equal to 0"
        self._unresolvedCodon = value
    
    @property
    def translationTable(self):
        return self._translationTable
    
    @translationTable.setter
    def translationTable(self, value):
        assert isinstance(value, int)
        assert value > 0, "translationTable must be greater than 0"
        self._translationTable = value
    
    def process(self):
        records = SeqIO.parse(open(self.fastaFile, 'r'), 'fasta')
        
        # Iterate over transcript records
        for record in records:
            # Declare output holding values that should reset for each transcript/record
            tempMProt = []
            tempMNucl = []
            tempAltProt = []
            tempAltNucl = []
            tempNoneProt = []
            tempNoneNucl = []
            
            # Iterate over strands
            for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                # Iterate over reading frames
                for frame in range(3):
                    length = 3 * ((len(record)-frame) // 3)
                    frameNuc = str(nuc[frame:frame+length])
                    frameProt = str(nuc[frame:frame+length].translate(table=self.translationTable))
                    
                    # Split protein/nucleotide into corresponding ORFs
                    ongoingLength = 0
                    splitNucleotide = []
                    splitProtein = []
                    frameProt = frameProt.split('*')
                    for i in range(len(frameProt)):
                        if len(frameProt) == 1 or i + 1 == len(frameProt):
                            splitProtein.append(frameProt[i])
                            splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])
                            ongoingLength += len(frameProt[i])*3
                        else:
                            splitProtein.append(frameProt[i] + '*')
                            splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                            ongoingLength += (len(frameProt[i]) + 1)*3
                    
                    # Fix unresolved regions
                    resolvedProt = []
                    resolvedNuc = []
                    indicesForDel = []
                    for i in range(len(splitProtein)):
                        if 'X' in splitProtein[i]:
                            posProt = []
                            for x in re.finditer(self.xRegex, splitProtein[i]):
                                if x.end() - x.start() > self.unresolvedCodon:
                                    posProt += [x.start(), x.end()]
                            if posProt == []:
                                continue
                            indicesForDel.insert(0, i)
                            
                            # Pull out resolved regions
                            resolvedProt.append(splitProtein[i][:posProt[0]])
                            resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                            for x in range(1, len(posProt)-1, 2):
                                start = posProt[x]
                                end = posProt[x+1]
                                resolvedProt.append(splitProtein[i][start:end])
                                resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                            resolvedProt.append(splitProtein[i][posProt[-1]:])
                            resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
                    
                    # Delete old entries and add resolved entries
                    for index in indicesForDel:
                        del splitProtein[index]
                        del splitNucleotide[index]
                    splitProtein += resolvedProt
                    splitNucleotide += resolvedNuc
                    
                    # Enter the main processing loop with our resolved regions
                    for i in range(len(splitProtein)):                              # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                        # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                        mPro = None
                        altPro = None
                        nonePro = None
                        codonIndex = None
                        noneCodonContingency = None
                        
                        # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                        if len(splitProtein[i]) < self.minProLen:            # Disregard sequences that won't meet the size requirement
                            continue
                        acceptedPro = str(splitProtein[i])
                        
                        # Alternative start coding      
                        nucSeqOfProt = splitNucleotide[i]               # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                        codons = re.findall('..?.?', nucSeqOfProt)          # Pulls out a list of codons from the nucleotide
                        for codon in codons:                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                            if codon == 'GTG' or codon == 'TTG':
                                codonIndex = codons.index(codon)    # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                break
                            elif codon == 'CTG':
                                if noneCodonContingency == None:    # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                    noneCodonContingency = codons.index(codon)
                        
                        # Get the three ORF versions from each region inbetween stop codons
                        if 'M' in str(acceptedPro):                 # Obtains a traditional methionine initiated ORF starting from the first methionine if there is one in the sequence
                            mPro = self.startCodon.search(str(acceptedPro)).groups()[0]  # Note that startCodon was declared at the start of this file     
                        
                        if codonIndex != None:                  # Gets the start position of the protein if we found a likely alternative start (aka a 'GTG' or 'TTG')
                            altPro = acceptedPro[codonIndex:]
                        elif noneCodonContingency != None:              # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                            altPro = acceptedPro[codonIndex:]
                        
                        nonePro = acceptedPro
                        
                        # Store if sequence lengths are within the desired range
                        if mPro != None:
                            if (len(mPro) >= self.minProLen) and (self.maxProLen == 0 or len(mPro) <= self.maxProLen):
                                tempMProt.append(mPro)
                                newStartPosition = acceptedPro.find(mPro)
                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                        if altPro != None:
                            if (len(altPro) >= self.minProLen) and (self.maxProLen == 0 or len(altPro) <= self.maxProLen):
                                tempAltProt.append(altPro)
                                newStartPosition = acceptedPro.find(altPro[1:]) - 1
                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                        if nonePro != None:
                            if (len(nonePro) >= self.minProLen) and (self.maxProLen == 0 or len(nonePro) <= self.maxProLen):
                                tempNoneProt.append(nonePro)
                                tempNoneNucl.append(str(nucSeqOfProt))
            
            # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
            if len(tempMProt + tempAltProt + tempNoneProt) >= 1:
                # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                ## Prot list
                for i in range(0, self.hitsToPull-len(tempMProt)):
                    tempMProt.append('-')
                for i in range(0, self.hitsToPull-len(tempAltProt)):
                    tempAltProt.append('-')
                for i in range(0, self.hitsToPull-len(tempNoneProt)):
                    tempNoneProt.append('-')
                ## Nucl list
                for i in range(0, self.hitsToPull-len(tempMNucl)):
                    tempMNucl.append('-')
                for i in range(0, self.hitsToPull-len(tempAltNucl)):
                    tempAltNucl.append('-')
                for i in range(0, self.hitsToPull-len(tempNoneNucl)):
                    tempNoneNucl.append('-')
                
                # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted)
                ## Prot
                tempSortedMProt = sorted(tempMProt, key=len)
                tempSortedAltProt = sorted(tempAltProt, key=len)
                tempSortedNoneProt = sorted(tempNoneProt, key=len)
                ## Nucl
                tempSortedMNucl = sorted(tempMNucl, key=len)
                tempSortedAltNucl = sorted(tempAltNucl, key=len)
                tempSortedNoneNucl = sorted(tempNoneNucl, key=len)
                
                # Run a final size comparison to choose the best ORF(s).
                tempOverallProt = []
                tempOverallNucl = []
                for i in range(0, self.hitsToPull):
                    # Accept a none protein if it's longer+cutoff than the alternatives
                    if len(tempSortedNoneProt[-1]) > len(tempSortedAltProt[-1]) + self.noCodonStringency and len(tempSortedNoneProt[-1]) > len(tempSortedMProt[-1]) + self.noCodonStringency:     # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                        tempOverallProt.append(tempSortedNoneProt[-1])
                        tempSortedNoneProt.pop()
                        
                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                        tempSortedNoneNucl.pop()
                    # Accept an alt protein if it's longer+cutoff than the methionine protein
                    elif len(tempSortedAltProt[-1]) > len(tempSortedMProt[-1]) + self.altCodonStringency:
                        tempOverallProt.append(tempSortedAltProt[-1])
                        tempSortedAltProt.pop()
                        
                        tempOverallNucl.append(tempSortedAltNucl[-1])
                        tempSortedAltNucl.pop()
                    # Accept an M protein if it isn't a blank
                    elif tempSortedMProt[-1] != "-":
                        tempOverallProt.append(tempSortedMProt[-1])                                                             # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                        tempSortedMProt.pop()
                        
                        tempOverallNucl.append(tempSortedMNucl[-1])                                                             # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                        tempSortedMNucl.pop()
                    # Default to a none protein if nothing else is available and it isn't blank
                    elif tempSortedNoneProt[-1] != "-":
                        tempOverallProt.append(tempSortedNoneProt[-1])
                        tempSortedNoneProt.pop()
                        
                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                        tempSortedNoneNucl.pop()
                
                # Yield results
                for i in range(0, self.hitsToPull):
                    if len(tempOverallProt) <= i:
                        break
                    if tempOverallProt[i] != "-":
                        yield record.id, tempOverallProt[i], tempOverallNucl[i]

class ORFPredictionProcess(BasicProcess):
    '''
    Handles ORF prediction in a separate thread.
    
    Parameters:
        mrnaFileIn -- a string indicating the location of the FASTA file containing mRNA sequences.
        cdsFileOut -- a string indicating the location to write the CDS sequences to.
        protFileOut -- a string indicating the location to write the protein sequences to.
    '''
    def task(self, mrnaFileIn, cdsFileOut, protFileOut, translationTable):
        # Establish the finder object
        orfFinder = ORF_Find(mrnaFileIn)
        orfFinder.hitsToPull = 1
        orfFinder.unresolvedCodon = 5 # allow some codons to be unresolved
        orfFinder.translationTable = translationTable
        
        # Produce the CDS and protein sequence files
        with open(cdsFileOut, "w") as cdsOut, open(protFileOut, "w") as protOut:
            for mrnaID, protSeq, cdsSeq in orfFinder.process():
                if protSeq == "-" or cdsSeq == "-": # '-' is a blank from ORF_Find
                    continue
                
                cdsOut.write(f">{mrnaID}\n{cdsSeq}\n")
                protOut.write(f">{mrnaID}\n{protSeq}\n")
        
        # Touch the OK files
        touch_ok(cdsFileOut)
        touch_ok(protFileOut)

# Define functions
def generate_sequence_length_index(fastaFile):
    '''
    Helper function to take a FASTA file and generate a dictionary of
    sequence lengths which is to be pickled at the same location as the FASTA.
    
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to index.
    '''
    seqLenDict = {}
    for title, seq in FastaParser(fastaFile):
        seqid = title.split(None, 1)[0]
        if seqid in seqLenDict:
            raise KeyError(f"'{seqid}' was found multiple times in '{fastaFile}'!")
        seqLenDict[seqid] = len(seq)
    
    indexFile = f"{fastaFile}.lengths.pkl"
    with open(indexFile, "wb") as fileOut:
        pickle.dump(seqLenDict, fileOut)

def remove_sequence_from_fasta(fastaFile, sequenceIDs, outputFile, force=False):
    '''
    Helper function to take a FASTA file and remove one or more sequences in the
    sequenceIDs list. The result is written to outputFile.
    
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to remove IDs from.
        sequenceIDs -- a list of strings indicating the sequence IDs to remove.
        outputFile -- a string indicating the location to write the modified FASTA to.
    '''
    # Check if output file exists
    if os.path.exists(outputFile):
        if force:
            os.unlink(outputFile)
        else:
            raise FileExistsError(f"remove_sequence_from_fasta() failed because '{outputFile}' exists!")
    
    # Iterate through the FASTA and drop any IDs encountered
    sequenceIDs = set(sequenceIDs)
    foundIDs = []
    wroteASequence = False
    with open(fastaFile, "r") as fileIn, open(outputFile, "w") as fileOut:
        records = SeqIO.parse(fileIn, "fasta")
        for record in records:
            seqID = record.id
            if seqID not in sequenceIDs:
                fileOut.write(record.format("fasta"))
                wroteASequence = True
            else:
                foundIDs.append(seqID)
    
    # Check if all IDs were found
    if len(sequenceIDs) != len(foundIDs):
        os.unlink(outputFile) # clean up flawed file
        missing = sequenceIDs - set(foundIDs)
        raise KeyError(("remove_sequence_from_fasta() failed because " +
                        f"the following IDs weren't found in '{fastaFile}': {missing}"))
    
    # Check if any sequences were written at all
    if not wroteASequence:
        os.unlink(outputFile)
        raise Exception(("remove_sequence_from_fasta() failed because no sequences " +
                         "remained after removing the provided sequence IDs!"))

def process_transcripts(locations, threads, translationTable=1):
    '''
    Will take the files within the 'transcripts' subdirectory of workingDirectory and
    produce sequence files from them where appropriate.
    
    Parameters:
        locations -- a Locations object with attributes for directory locations.
        threads -- an integer indicating the number of threads to use for ORF prediction.
        translationTable -- an integer indicating the translation table to use for ORF prediction;
                            defaults to 1 which is the standard table.
    '''
    # Locate all transcript files / triplets
    needsSymlink = []
    needsPrediction = []
    for file in os.listdir(locations.txDir):
        if file.endswith(".mrna"):
            if not file.startswith("transcriptome"):
                raise ValueError(f"'{file}' in '{locations.txDir}' does not begin with 'transcriptome' as expected")
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".mrna")[0]
            suffixNum = filePrefix.split("transcriptome")[1]
            if not suffixNum.isdigit():
                raise ValueError(f"'{file}' in '{locations.txDir}' does not have a number suffix as expected")
            
            # Symbolic link the mRNA file
            mrnaFile = os.path.join(locations.txDir, file)
            mrnaLink = os.path.join(locations.sequencesDir, file)
            if os.path.exists(mrnaLink):
                handle_symlink_change(mrnaLink, mrnaFile)
            else:
                os.symlink(mrnaFile, mrnaLink)
                touch_ok(mrnaLink)
            
            # Build up the triplet of files
            thisTriplet = [mrnaFile, None, None]
            cdsFile = os.path.join(locations.txDir, f"transcriptome{suffixNum}.cds")
            aaFile = os.path.join(locations.txDir, f"transcriptome{suffixNum}.aa")
            if os.path.exists(cdsFile) and os.path.exists(aaFile):
                thisTriplet[1] = cdsFile
                thisTriplet[2] = aaFile
            
            # Store the triplet according to their needs
            if None in thisTriplet:
                needsPrediction.append([thisTriplet, suffixNum])
            else:
                needsSymlink.append([thisTriplet, suffixNum])
    
    # Symlink files that don't need ORF prediction
    for (mrnaFile, cdsFile, protFile), suffixNum in needsSymlink:
        # Derive output file names
        mrnaLink = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.mrna")
        cdsLink = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.cds")
        protLink = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.aa")
        
        # Symlink the files
        for origFile, linkFile in zip([cdsFile, protFile], [cdsLink, protLink]):
            if os.path.exists(linkFile):
                handle_symlink_change(linkFile, origFile)
            else:
                os.symlink(origFile, linkFile)
                touch_ok(linkFile)
    
    # Narrow down files needing prediction to those not previously predicted
    filteredPrediction = []
    for triplet, suffixNum in needsPrediction:
        # Derive output file names
        cdsFileName = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.cds")
        protFileName = os.path.join(locations.sequencesDir, f"transcriptome{suffixNum}.aa")
        newTriplet = [cdsFileName, protFileName]
        
        # Store files that need ORF prediction
        if not all([ os.path.exists(x) for x in newTriplet ]) and not all([ os.path.exists(x + ".ok") for x in newTriplet ]):
            filteredPrediction.append([triplet[0], newTriplet])
    
    # Predict ORFs for files needing it
    if len(filteredPrediction) > 0:
        print(f"# Predicting ORFs from transcriptome file(s)...")
        
        for i in range(0, len(filteredPrediction), threads): # only process n (threads) files at a time
            processing = []
            for x in range(threads): # begin processing n files
                if i+x < len(filteredPrediction): # parent loop may excess if n > the number of files needing indexing
                    mrnaFile, (cdsFile, protFile) = filteredPrediction[i+x]
                    
                    predictionWorkerThread = ORFPredictionProcess(mrnaFile, cdsFile, protFile, translationTable)
                    predictionWorkerThread.start()
                    processing.append(predictionWorkerThread)
            
            # Gather results
            for predictionWorkerThread in processing:
                predictionWorkerThread.join()
                predictionWorkerThread.check_errors()
