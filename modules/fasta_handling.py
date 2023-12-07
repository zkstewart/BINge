import os, sys, re
from Bio import SeqIO
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_SeqIO # expose this to any callers

class AnnotationExtractor:
    '''
    Class to encapsulate the logic of taking a GFF3, its associated FASTA, and generating
    sequences for mRNA features. It is designed for use with BINge as a (hopefully) light
    weight actor.
    '''
    def __init__(self, gff3File, fastaFile):
        self.gff3File = gff3File
        self._fastaFile = fastaFile
        
        self.fasta = SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))
        
        self.parentFeatureType = "gene"
        self.childFeatureType = "mRNA"
    
    def _gff3_iterator(self):
        '''
        Provides a simple iterator for a GFF3 file that AnnotationExtractor is built around.
        Will only parse features where the parent and child match the object attributes.
        
        This parser is intended to be quick and dumb - which is a positive for GFF3s which are
        frequently poorly formatted. But if a GFF3 is _so poorly_ formatted that it isn't even
        ordered as expected (gene -> mRNA -> CDS/exon, then the next gene -> ...) this function
        will fail.
        '''
        idRegex = re.compile(r"ID=(.+?)($|;|\n)")
        details = [None, None, None] # mrnaID, contig, strand
        exon = []
        cds = []
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
                isMRNA = featureType == "mRNA"
                
                # Check if we should yield a feature
                if isMRNA and details[0] != None:
                    if len(exon) > 0 and len(cds) > 0: # if this fails it might be a pseudogene
                        yield details, exon, cds
                
                # Check if we should build a new feature
                if isMRNA:
                    details = [idRegex.search(sl[8]).groups()[0], contig, strand]
                    exon = []
                    cds = []
                
                # Build an ongoing feature
                if featureType == "exon":
                    exon.append([int(start), int(end)])
                elif featureType == "CDS":
                    cds.append([int(start), int(end), int(frame)])
        
        # Yield the last feature in the GFF3
        if details != [None, None, None]:
            yield details, exon, cds
        else:
            raise Exception(f"AnnotationExtractor error: '{self.gff3File}' does not appear to be a valid GFF3 file!")
    
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
            ".assemble_sequence() received coordinates that aren't in a list of lists with 2 or 3 entries each!"
        
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
        for details, exon, cds in self._gff3_iterator():
            mrnaID, contig, strand = details
            
            # Create sequence by piecing together exon / CDS bits
            contigSequence = str(self.fasta[contig].seq)
            
            exonSeq = AnnotationExtractor.assemble_sequence(exon, strand, contigSequence)
            cdsSeq = AnnotationExtractor.assemble_sequence(cds, strand, contigSequence)
            
            # Yield products
            yield mrnaID, exonSeq, cdsSeq

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
