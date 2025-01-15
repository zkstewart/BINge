import os, sys, re, pickle
from Bio import SeqIO
from pyfaidx import Fasta
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .thread_workers import BasicProcess
from .validation import handle_symlink_change, touch_ok

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_ORF

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
        orfFinder = ZS_ORF.ORF_Find(mrnaFileIn)
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
