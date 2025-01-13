import os, sys, re, pickle
from Bio import SeqIO
from pyfaidx import Fasta
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .thread_workers import BasicProcess
from .validation import handle_symlink_change, touch_ok

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_SeqIO, ZS_ORF

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
        
        self._parse_fastas()
    
    def _parse_fastas(self):
        for fastaFile in self.fastaFiles:
            self.records.append(Fasta(fastaFile))
    
    def __getitem__(self, key):
        for records in self.records:
            try:
                return records[key]
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

class AnnotationExtractor:
    '''
    Class to encapsulate the logic of taking a GFF3, its associated FASTA, and generating
    sequences for mRNA features. It is designed for use with BINge as a (hopefully) light
    weight actor.
    '''
    def __init__(self, gff3File, fastaFile, isMicrobial=False):
        self.gff3File = gff3File
        self.isMicrobial = isMicrobial
        self._fastaFile = fastaFile
        
        self.fasta = SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))
    
    def _gff3_iterator(self):
        '''
        Provides a simple iterator for a GFF3 file that AnnotationExtractor is built around.
        
        It is intended to be quick and dumb - which is a positive for GFF3s which are
        frequently poorly formatted. But if a GFF3 is _so poorly_ formatted that it isn't even
        ordered as expected (gene -> mRNA -> CDS/exon, then the next gene -> ...) this function
        will fail.
        
        These failures cannot be detected so the user must be vigilant for them.
        '''
        # Setup data structures
        idRegex = re.compile(r"ID=(.+?)($|;|\n)")
        parentRegex = re.compile(r"Parent=(.+?)($|;|\n)")
        foundSomething = False
        gff3Dict = {}
        
        # Configuration for microbial or normal GFF3
        if self.isMicrobial:
            subfeature = ["gene"]
        else:
            subfeature = ["mRNA", "transcript"]
        
        # Iterate through file and build up a dictionary of features
        with open(self.gff3File, "r") as fileIn:
            for line in fileIn:
                # Parse and skip irrelevant lines
                if line.startswith("#"):
                    continue
                sl = line.rstrip("\t\r\n ").split("\t")
                if len(sl) != 9: # gff3 lines are always 9 long
                    continue
                
                # Extract line details
                contig, source, featureType, start, end, \
                    score, strand, frame, attributes = sl
                isParent = featureType in subfeature
                
                # Index parent features (which may be genes or mRNAs)
                if isParent:
                    mrnaID = idRegex.search(sl[8]).groups()[0]
                    gff3Dict[mrnaID] = [contig, strand, [], []] # exon, CDS lists
                
                # Build an ongoing feature
                featureIndex = 2 if featureType == "exon" else 3 if featureType == "CDS" else None
                if featureIndex != None:
                    parentIDs = parentRegex.search(sl[8]).groups()[0].split(",")
                    for pid in parentIDs:
                        if pid in gff3Dict:
                            gff3Dict[pid][featureIndex].append([int(start), int(end), frame])
                            foundSomething = True
        
        # Yield each feature in the GFF3
        if foundSomething == False:
            raise Exception(f"AnnotationExtractor error: '{self.gff3File}' does not appear to be a valid GFF3 file!")
        else:
            for mrnaID, (contig, strand, exon, cds) in gff3Dict.items():
                if not self.isMicrobial:
                    yield mrnaID, contig, strand, exon, cds
                else:
                    yield mrnaID, contig, strand, cds, cds # exon is equivalent to CDS in microbe GFF3s
    
    @staticmethod
    def assemble_sequence(coordsList, strand, contigSequence):
        '''
        Function to receive a list of exon or CDS coordinates and assemble those into a nucleotide sequence.
        
        Parameters:
            coordsList -- a list of lists containing [start, end] value if exon, or [start, end, frame] if CDS
            strand -- either "+" or "-" for positive and negative strands, respectively
            contigSequence -- a string of the sequence to extract bits out of.
        '''
        assert all([ len(x) == 2 for x in coordsList ]) or all([ len(x) == 3 for x in coordsList ]), \
            "assemble_sequence() received coordinates that aren't in a list of lists with 2 or 3 entries each!"
        
        # Sort coordsList appropriately
        coordsList.sort(key = lambda x: (int(x[0]), int(x[1])))
        
        # Assemble the sequence now
        sequence = ""
        for i, value in enumerate(coordsList):
            # Extract details from coords value
            if len(value) == 2:
                start, end = value
                frame = 0
            else:
                start, end, frame = value
            
            # Extend sequence
            sequenceBit = contigSequence[start-1:end] # 1-based correction to start to make it 0-based
            sequence += sequenceBit
            
            # Note the frame if necessary
            if (i == 0 and strand == "+") or (i == len(coordsList)-1 and strand == "-"):
                try:
                    startingFrame = int(frame)
                except:
                    startingFrame = 0
        
        # Reverse complement if necessary
        if strand == "-":
            sequence = ZS_SeqIO.FastASeq.get_reverse_complement(self=None, staticSeq=sequence)
        
        # Trim the sequence to the correct frame
        sequence = sequence[startingFrame:]
        
        return sequence
    
    def iter_sequences(self):
        '''
        Iterates through the GFF3 file and yields CDS and exon sequences associated with each mRNA feature.
        
        Yields:
            mrnaID -- a string of the mRNA sequence's identifier.
            exonSeq -- a string representing the exon/transcript sequence.
            cdsSeq -- a string representing the CDS.
        '''
        warnedOnce1 = False
        warnedOnce2 = False
        for mrnaID, contig, strand, exon, cds in self._gff3_iterator():
            # Skip if the CDS or exon lists are empty
            "This can happen with 'transcript' biotypes that don't have CDS features"
            if len(exon) == 0 or len(cds) == 0:
                if not warnedOnce1:
                    print(f"WARNING: '{mrnaID}' lacks exon and/or CDS features; " + 
                          f"similar warnings for file '{self.gff3File}' will be suppressed.")
                    warnedOnce1 = True
                continue
            
            # Create sequence by piecing together exon / CDS bits
            contigSequence = str(self.fasta[contig].seq)
            
            exonSeq = AnnotationExtractor.assemble_sequence(exon, strand, contigSequence)
            cdsSeq = AnnotationExtractor.assemble_sequence(cds, strand, contigSequence)
            protSeq = ZS_SeqIO.FastASeq.dna_to_protein(cdsSeq)
            
            # Warn if the protein sequence contains a stop codon
            if "*" in protSeq.strip("*"):
                if not warnedOnce2:
                    print(f"WARNING: protein sequence for '{mrnaID}' contains an internal stop codon; " + 
                          f"similar warnings for file '{self.gff3File}' will be suppressed.")
                    warnedOnce2 = True
            
            # Yield products
            yield mrnaID, exonSeq, cdsSeq, protSeq

class ORFPredictionProcess(BasicProcess):
    '''
    Handles ORF prediction in a separate thread.
    
    Parameters:
        mrnaFileIn -- a string indicating the location of the FASTA file containing mRNA sequences.
        cdsFileOut -- a string indicating the location to write the CDS sequences to.
        protFileOut -- a string indicating the location to write the protein sequences to.
    '''
    def task(self, mrnaFileIn, cdsFileOut, protFileOut):
        # Establish the finder object
        orfFinder = ZS_ORF.ORF_Find(mrnaFileIn)
        orfFinder.hitsToPull = 1
        orfFinder.minProLen = 1 # we need a hit for _every_ sequence to prevent errors downstream
        orfFinder.unresolvedCodon = 5 # allow some codons to be unresolved
        
        # Produce the CDS and protein sequence files
        with open(cdsFileOut, "w") as cdsOut, open(protFileOut, "w") as protOut:
            for mrnaID, protSeq, cdsSeq in orfFinder.process():
                if protSeqs == "-" or cdsSeq == "-": # '-' is a blank from ORF_Find
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

def process_transcripts(locations, threads):
    '''
    Will take the files within the 'transcripts' subdirectory of workingDirectory and
    produce sequence files from them where appropriate.
    
    Parameters:
        locations -- a Locations object with attributes for directory locations.
        threads -- an integer indicating the number of threads to use for ORF prediction.
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
                    
                    predictionWorkerThread = ORFPredictionProcess(mrnaFile, cdsFile, protFile)
                    predictionWorkerThread.start()
                    processing.append(predictionWorkerThread)
            
            # Gather results
            for predictionWorkerThread in processing:
                predictionWorkerThread.join()
                predictionWorkerThread.check_errors()
