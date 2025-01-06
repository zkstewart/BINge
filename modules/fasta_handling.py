import os, sys, re, pickle
from Bio import SeqIO
from pyfaidx import Fasta
from Bio.SeqIO.FastaIO import SimpleFastaParser

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_SeqIO, ZS_ORF # expose this to any callers

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
        Will only parse features where the parent and child match the object attributes.
        
        This parser is intended to be quick and dumb - which is a positive for GFF3s which are
        frequently poorly formatted. But if a GFF3 is _so poorly_ formatted that it isn't even
        ordered as expected (gene -> mRNA -> CDS/exon, then the next gene -> ...) this function
        will fail.
        
        These failures will be detected and the program will end to prevent erroneous behaviour.
        '''
        # Setup data structures
        idRegex = re.compile(r"ID=(.+?)($|;|\n)")
        parentRegex = re.compile(r"Parent=(.+?)($|;|\n)")
        details = [None, None, None] # mrnaID, contig, strand
        exon = []
        cds = []
        foundSomething = False
        
        # Configuration for microbial or normal GFF3
        if self.isMicrobial:
            subfeature = ["gene"]
        else:
            subfeature = ["mRNA", "transcript"]
        
        # Iterate through file now
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
                    score, strand, frame, attributes \
                    = sl
                isMRNA = featureType in subfeature
                isBreakpoint = featureType == "gene"
                
                # Check if we should yield a feature
                if (isMRNA or isBreakpoint) and details[0] != None:
                    if not self.isMicrobial:
                        if len(exon) > 0 and len(cds) > 0: # if this fails it might be a pseudogene
                            foundSomething = True
                            yield details, exon, cds
                    else:
                        if len(cds) > 0: # if this fails it might be a pseudogene
                            foundSomething = True
                            yield details, cds, cds # CDS and exon are equivalent in microbe GFF3s
                    
                    # Zero the details again as a safeguard for if we hit a breakpoint and not isMRNA
                    details = [None, None, None]
                    exon = []
                    cds = []
                
                # Check if we should build a new feature
                if isMRNA:
                    details = [idRegex.search(sl[8]).groups()[0], contig, strand]
                    exon = []
                    cds = []
                
                # Build an ongoing feature
                if featureType == "exon":
                    if details[0] != None:
                        assert parentRegex.search(sl[8]).groups()[0] == details[0], \
                            "AnnotationExtractor error: exon parent doesn't match mRNA ID; file is not ordered!!"
                    exon.append([int(start), int(end)])
                elif featureType == "CDS":
                    if details[0] != None:
                        assert parentRegex.search(sl[8]).groups()[0] == details[0], \
                            "AnnotationExtractor error: CDS parent doesn't match mRNA ID; file is not ordered!!"
                    cds.append([int(start), int(end), int(frame)])
        
        # Yield the last feature in the GFF3
        if foundSomething == False:
            raise Exception(f"AnnotationExtractor error: '{self.gff3File}' does not appear to be a valid GFF3 file!")
        if details[0] != None:
            if not self.isMicrobial:
                yield details, exon, cds
            else:
                yield details, cds, cds
    
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
        for value in coordsList:
            if len(value) == 2:
                start, end = value
            else:
                start, end, _ = value # don't need the frame value for this
            
            sequenceBit = contigSequence[start-1:end] # 1-based correction to start to make it 0-based
            sequence += sequenceBit
        
        # Reverse complement if necessary
        if strand == "-":
            sequence = ZS_SeqIO.FastASeq.get_reverse_complement(self=None, staticSeq=sequence)
        
        return sequence
    
    def iter_sequences(self):
        '''
        Iterates through the GFF3 file and yields CDS and exon sequences associated with each mRNA feature.
        
        Yields:
            mrnaID -- a string of the mRNA sequence's identifier.
            exonSeq -- a string representing the exon/transcript sequence.
            cdsSeq -- a string representing the CDS.
        '''
        warnedOnce = False
        for details, exon, cds in self._gff3_iterator():
            mrnaID, contig, strand = details
            
            # Create sequence by piecing together exon / CDS bits
            contigSequence = str(self.fasta[contig].seq)
            
            exonSeq = AnnotationExtractor.assemble_sequence(exon, strand, contigSequence)
            cdsSeq = AnnotationExtractor.assemble_sequence(cds, strand, contigSequence)
            protSeq = ZS_SeqIO.FastASeq.dna_to_protein(cdsSeq)
            
            # Warn if the protein sequence contains a stop codon
            if "*" in protSeq.strip("*"):
                if not warnedOnce:
                    print(f"WARNING: protein sequence for '{mrnaID}' contains an internal stop codon; " + 
                          f"further warnings for file '{self.gff3File}' will be suppressed.")
                    warnedOnce = True
            
            # Yield products
            yield mrnaID, exonSeq, cdsSeq, protSeq

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

def load_sequence_length_index(indexFile):
    '''
    Load in the pickled result of generate_sequence_length_index().
    
    Parameters:
        indexFile -- a string indicating the location of index generated from a FASTA file.
    '''
    if os.path.exists(indexFile):
        with open(indexFile, "rb") as fileIn:
            seqLenDict = pickle.load(fileIn)
        return seqLenDict
    else:
        raise FileNotFoundError(("load_sequence_length_index() failed because " + 
                                 f"'{indexFile}' doesn't exist!"))

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
