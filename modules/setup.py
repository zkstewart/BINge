import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from validation import validate_fasta, handle_symlink_change, touch_ok
from fasta_handling import generate_sequence_length_index, txome_to_orfs
from gmap_handling import build_index
from gff3tofasta import gff3_to_fasta

def linker(directory, prefix, suffix, file):
    symlink = os.path.join(directory, f"{prefix}.{suffix}" + (".gz" if file.endswith(".gz") else ""))
    if os.path.exists(symlink):
        handle_symlink_change(symlink, file)
    else:
        os.symlink(file, symlink)
    return symlink

def gff3_extractor(extractionDirectory, fasta, gff3, prefix, isMicrobial, translationTable):
    # Derive names for output sequences
    mrna = os.path.join(extractionDirectory, f"{prefix}.mrna")
    cds = os.path.join(extractionDirectory, f"{prefix}.cds")
    aa = os.path.join(extractionDirectory, f"{prefix}.aa")
    outputFileNames = [mrna, cds, aa]
    
    # Check for file existence
    filesExist = all([ os.path.exists(x) for x in outputFileNames ])
    flagsExist = all([ os.path.exists(x + ".ok") for x in outputFileNames ])
    
    # Extract annotations if applicable
    if not (filesExist and flagsExist):
        print(f"# Extracting sequences from '{prefix}'")
        gff3_to_fasta(gff3, fasta, mrna, cds, aa, isMicrobial, translationTable)
    
    return mrna, cds, aa

def orf_extractor(extractionDirectory, inputFile, prefix, translationTable):
    # Derive names for output sequences
    mrna = os.path.join(extractionDirectory, f"{prefix}.mrna")
    cds = os.path.join(extractionDirectory, f"{prefix}.cds")
    aa = os.path.join(extractionDirectory, f"{prefix}.aa")
    outputFileNames = [mrna, cds, aa]
    
    # Check for file existence
    filesExist = all([ os.path.exists(x) for x in outputFileNames ])
    flagsExist = all([ os.path.exists(x + ".ok") for x in outputFileNames ])
    
    # Extract ORFs if applicable
    if not (filesExist and flagsExist):
        print(f"# Extracting sequences from '{prefix}'")
        txome_to_orfs(inputFile, cds, aa, translationTable)
        
        # Symlink the original input file
        linker(extractionDirectory, prefix, "mrna", inputFile)
        touch_ok(mrna) # txome_to_orfs() touches ok files for cds and aa
    
    return cds, aa

class TargetGenome:
    '''
    Class assumes that input files are pre-validated to exist.
    '''
    def __init__(self, index, locations, fasta, gff3=None):
        self.prefix = f"genome{index}"
        self.directory = locations.genomesDir
        self.extractionDirectory = locations.genomesDir
        self.gff3 = gff3 # enter this value first
        self.fasta = fasta # so this can raise a better error msg if needed
        
        self.thisType = "target"
        self.mrna = None
        self.cds = None
        self.aa = None
    
    @property
    def gmap_db(self):
        return os.path.join(self.directory, f"{self.prefix}.gmap")
    
    @property
    def fasta(self):
        return self._fasta # this is the symlink
    
    @fasta.setter
    def fasta(self, value):
        isFASTA = validate_fasta(value)
        if not isFASTA:
            if self.gff3 != None:
                raise ValueError(f"-i value '{value}' is not a FASTA file; make sure you order " +
                                 "input as gff3,fasta")
            else:
                raise ValueError(f"-i value '{value}' is not a FASTA file")
        
        symlink = linker(self.directory, self.prefix, "fasta", os.path.abspath(value))
        self._fasta = symlink
    
    @property
    def gff3(self):
        return self._gff3 # this is the symlink
    
    @gff3.setter
    def gff3(self, value):
        if value == None:
            self._gff3 = None
        else:
            isFASTA = validate_fasta(value)
            if isFASTA:
                raise ValueError(f"-i value '{value}' is not a GFF3 file; make sure you order " +
                                 "input as gff3,fasta")
            
            symlink = linker(self.directory, self.prefix, "gff3", os.path.abspath(value))
            self._gff3 = symlink
    
    def extract_sequences(self, isMicrobial, translationTable):
        if self.gff3 != None:
            self.mrna, self.cds, self.aa = gff3_extractor(self.extractionDirectory, self.fasta,
                                                          self.gff3, self.prefix,
                                                          isMicrobial, translationTable)
        return self # needed to return this modified object back through the ProcessPoolExecutor
    
    def length_index(self):
        lengthFile = os.path.join(self.directory, f"{self.prefix}.lengths.pkl")
        if not (os.path.exists(lengthFile) and os.path.exists(f"{lengthFile}.ok")):
            generate_sequence_length_index(self.fasta, lengthFile)
            touch_ok(lengthFile)
    
    def gmap_index(self, gmapDir):
        #indexDir = os.path.join(self.directory, f"{self.prefix}.gmap")
        if not (os.path.exists(self.gmap_db) and os.path.exists(f"{self.gmap_db}.ok")):
            print(f"# Building GMAP index for '{self.prefix}'")
            build_index(os.path.join(gmapDir, "gmap_build"), # this location should be validated beforehand
                        self.fasta, self.gmap_db)
            touch_ok(self.gmap_db)
    
    def __repr__(self):
        return (f"<TargetGenome;prefix={self.prefix};fasta={self.fasta};gff3={self.gff3};" + 
               f"mrna={self.mrna};cds={self.cds};aa={self.aa}>")

class AnnotatedGenome:
    def __init__(self, index, locations, fasta, gff3):
        self.prefix = f"annotations{index}"
        self.directory = locations.gff3Dir
        self.extractionDirectory = locations.sequencesDir
        self.gff3 = gff3 # order of entering does not matter
        self.fasta = fasta # here as it did in TargetGenome
        
        self.thisType = "annotated"
        self.mrna = None
        self.cds = None
        self.aa = None
    
    @property
    def fasta(self):
        return self._fasta # this is the symlink
    
    @fasta.setter
    def fasta(self, value):
        isFASTA = validate_fasta(value)
        if not isFASTA:
            raise ValueError(f"--ig value '{value}' is not a FASTA file; make sure you order " +
                             "input as gff3,fasta")
        
        symlink = linker(self.directory, self.prefix, "fasta", os.path.abspath(value))
        self._fasta = symlink
    
    @property
    def gff3(self):
        return self._gff3 # this is the symlink
    
    @gff3.setter
    def gff3(self, value):
        if value == None:
            self._gff3 = None
        else:
            isFASTA = validate_fasta(value)
            if isFASTA:
                raise ValueError(f"--ig value '{value}' is not a GFF3 file; make sure you order " +
                                 "input as gff3,fasta")
            
            symlink = linker(self.directory, self.prefix, "gff3", os.path.abspath(value))
            self._gff3 = symlink
    
    def extract_sequences(self, isMicrobial, translationTable):
        self.mrna, self.cds, self.aa = gff3_extractor(self.extractionDirectory, self.fasta,
                                                      self.gff3, self.prefix,
                                                      isMicrobial, translationTable)
        return self # needed to return this modified object back through the ProcessPoolExecutor
    
    def __repr__(self):
        return f"<AnnotatedGenome;prefix={self.prefix};fasta={self.fasta};gff3={self.gff3}>"

class Transcriptome:
    def __init__(self, index, locations, mrnaFasta, cdsFasta=None, protFasta=None):
        self.prefix = f"transcriptome{index}"
        self.directory = locations.txDir
        self.extractionDirectory = locations.sequencesDir
        
        self.thisType = "transcriptome"
        self.mrna = mrnaFasta
        self.cds = cdsFasta
        self.aa = protFasta
    
    @property
    def mrna(self):
        '''
        At first this will be a symlink in self.directory. After running self.extract_sequences()
        it will instead be a symlink in self.extractionDirectory.
        '''
        return self._mrna
    
    @mrna.setter
    def mrna(self, value):
        isFASTA = validate_fasta(value)
        if not isFASTA:
            raise ValueError(f"--ix value '{value}' is not a FASTA file")
        
        symlink = linker(self.directory, self.prefix, "mrna", os.path.abspath(value))
        self._mrna = symlink
    
    @property
    def cds(self):
        return self._cds # this is the symlink
    
    @cds.setter
    def cds(self, value):
        if value == None:
            self._cds = None
        else:
            isFASTA = validate_fasta(value)
            if not isFASTA:
                raise ValueError(f"--ix value '{value}' is not a FASTA file")
            
            symlink = linker(self.directory, self.prefix, "cds", os.path.abspath(value))
            self._cds = symlink
    
    @property
    def aa(self):
        return self._aa # this is the symlink
    
    @aa.setter
    def aa(self, value):
        if value == None:
            self._aa = None
        else:
            isFASTA = validate_fasta(value)
            if not isFASTA:
                raise ValueError(f"--ix value '{value}' is not a FASTA file")
            
            symlink = linker(self.directory, self.prefix, "aa", os.path.abspath(value))
            self._aa = symlink
    
    def extract_sequences(self, translationTable):
        '''
        Unlike the other extract_sequences() functions, this one will not set self.mrna
        since doing so would 1) be unnecessary, as the file does not change. And 2) cause
        weirdness with the automatic symlinking built into the self.mrna property. Instead,
        orf_extractor() is responsible for setting the new symlink.
        '''
        if self.cds == None or self.aa == None:
            self.cds, self.aa = orf_extractor(self.extractionDirectory, self.mrna, ## TBD: Symlinks shouldn't be set by orf_extractor
                                              self.prefix, translationTable)
        return self # needed to return this modified object back through the ProcessPoolExecutor
    
    def __repr__(self):
        return f"<Transcriptome;prefix={self.prefix};mrna={self.mrna};cds={self.cds};aa={self.aa}>"
