import os, sys, re
from .validation import touch_ok
from .thread_workers import BasicProcess

# Define classes
class AnnotationExtractor:
    '''
    Class to encapsulate the logic of taking a GFF3, its associated FASTA, and generating
    sequences for mRNA features. It is designed for use with BINge as a (hopefully) light
    weight actor.
    '''
    def __init__(self, gff3File, fastaFile, isMicrobial=False, translationTable=1):
        self.gff3File = gff3File
        self.isMicrobial = isMicrobial
        self.translationTable = translationTable
        self._fastaFile = fastaFile
        
        self.fasta = SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))
    
    @property
    def translationTable(self):
        return self._translationTable
    
    @translationTable.setter
    def translationTable(self, value):
        assert isinstance(value, int)
        assert value > 0, "translationTable must be greater than 0"
        self._translationTable = value
    
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
            sequence = get_reverse_complement(sequence)
        
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
        for mrnaID, contig, strand, exon, cds in gff3_iterator(self.gff3File, self.isMicrobial):
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
            protSeq = dna_to_protein(cdsSeq, self.translationTable)
            
            # Warn if the protein sequence contains a stop codon
            if "*" in protSeq.strip("*"):
                if not warnedOnce2:
                    print(f"WARNING: protein sequence for '{mrnaID}' contains an internal stop codon; " + 
                          f"similar warnings for file '{self.gff3File}' will be suppressed.")
                    warnedOnce2 = True
            
            # Yield products
            yield mrnaID, exonSeq, cdsSeq, protSeq

# Multithreaded functions and classes
class GFF3ExtractionProcess(BasicProcess):
    '''
    Handles sequence extraction from a GFF3 file and genome FASTA file.
    
    Parameters:
        gff3FileIn -- a string indicating the location of a GFF3 formatted file containing annotations
                      corresponding to a genome.
        fastaFileIn -- a string indicating the location of a FASTA formatted file containing the genome
                       sequences that pair with the GFF3 file.
        mrnaFileOut -- a string indicating the location to write the mRNA sequences to.
        cdsFileOut -- a string indicating the location to write the CDS sequences to.
        protFileOut -- a string indicating the location to write the protein sequences to.
        isMicrobial -- a boolean indicating whether the organism is a microbe and hence GFF3
                       has gene -> CDS features, rather than gene -> mRNA -> CDS/exon.
    '''
    def task(self, gff3FileIn, fastaFileIn, mrnaFileOut, cdsFileOut, protFileOut, isMicrobial, translationTable):
        # Generate sequences
        with open(mrnaFileOut, "w") as mrnaOut, open(cdsFileOut, "w") as cdsOut, open(protFileOut, "w") as protOut:
            seqGenerator = AnnotationExtractor(gff3FileIn, fastaFileIn, isMicrobial, translationTable)
            for mrnaID, exonSeq, cdsSeq, protSeq in seqGenerator.iter_sequences():
                mrnaOut.write(f">{mrnaID}\n{exonSeq}\n")
                cdsOut.write(f">{mrnaID}\n{cdsSeq}\n")
                protOut.write(f">{mrnaID}\n{protSeq}\n")
        
        # Touch the .ok files to indicate completion
        touch_ok(mrnaFileOut)
        touch_ok(cdsFileOut)
        touch_ok(protFileOut)

# Other functions
def _build_data_to_yield(thisFeature):
    # Parse out variables from mRNA line
    contig, _, _, start, end, \
        _, strand, _, attributes \
        = thisFeature[1]
    
    # Start building a data structure to yield
    dataDict = {"contig": contig, "start": int(start), "end": int(end), "strand": strand}
    
    # Get the GMAP alignment quality attributes
    for attribute in attributes.split(";"):
        key, value = attribute.split("=")
        if key in ["identity", "coverage"]:
            dataDict[key] = float(value)
        elif key == "indels":
            dataDict[key] = int(value)
        elif key == "Name":
            dataDict["Name"] = value
    
    # Get the exons
    dataDict["exons"] = []
    for _sl in thisFeature[2:]:
        contig, _, featureType, start, end, \
            _, strand, _, attributes \
            = _sl
        if featureType == "exon":
            dataDict["exons"].append([int(start), int(end)])
    dataDict["exons"] = sorted(dataDict["exons"])
    
    return dataDict

def iterate_gmap_gff3(gff3File):
    '''
    Provides a simple iterator for a GMAP GFF3 file which yields only the
    relevant details for BINge. Assumes the file is sorted as is the case
    with GMAP output.
    
    Parameters:
        gff3File -- a string indicating the location of a GMAP GFF3 formatted
                    file that is sorted.
    '''
    thisAlignment = []
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\t\r\n ").split("\t")
            
            # Skip comment / irrelevant lines
            if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                continue
            
            # Handle content lines
            else:
                featureType = sl[2]
                
                # Yield a completed gene data dictionary
                if featureType == "gene" and thisAlignment != []:
                    dataDict = _build_data_to_yield(thisAlignment)
                    
                    # Reset our feature storage, and yield result now
                    thisAlignment = []
                    yield dataDict
                
                # Build an ongoing feature
                thisAlignment.append(sl)
    
    # Yield the last feature in the GFF3
    if thisAlignment != []:
        dataDict = _build_data_to_yield(thisAlignment)
    else:
        raise Exception(f"'{gff3File}' does not appear to be a valid GFF3 file!")
    yield dataDict

def gff3_iterator(gff3File, isMicrobial=False):
    '''
    Provides a simple iterator for a GFF3 file.
    
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
    if isMicrobial:
        subfeature = ["gene"]
    else:
        subfeature = ["mRNA", "transcript"]
    
    # Iterate through file and build up a dictionary of features
    with open(gff3File, "r") as fileIn:
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
        raise Exception(f"'{gff3File}' does not appear to be a valid GFF3 file!")
    else:
        for mrnaID, (contig, strand, exon, cds) in gff3Dict.items():
            if not isMicrobial:
                yield mrnaID, contig, strand, exon, cds
            else:
                yield mrnaID, contig, strand, cds, cds # exon is equivalent to CDS in microbe GFF3s

def extract_annotations_from_gff3(inputDir, outputDir, outputPrefix, threads, isMicrobial, translationTable):
    '''
    Will take the files within the inputDir location and produce sequence files from them
    which are written to the outputDir location.
    
    Parameters:
        inputDir -- a string indicating the location of a directory containing GFF3/FASTA files.
        outputDir -- a string indicating the location to write the output sequence files to.
        outputPrefix -- a string indicating the prefix to use for the output sequence files.
        threads -- an integer indicating the number of threads to use for processing.
        isMicrobial -- a boolean indicating whether the organism is a microbe and hence GFF3
                       has gene -> CDS features, rather than gene -> mRNA -> CDS/exon.
        translationTable -- an integer indicating the translation table to use for protein sequences;
                            should be an integer corresponding to the NCBI translation table numbering.
    '''
    suffixRegex = re.compile(r"(\d+)\.gff3$")
    
    # Locate all GFF3/genome pairs
    filePairs = []
    for file in os.listdir(inputDir):
        if file.endswith(".gff3"):
            # Extract suffix number from file name
            suffixNum = suffixRegex.search(file).group(1)
            if not suffixNum.isdigit():
                raise ValueError(f"'{file}' in '{inputDir}' does not have a number suffix as expected")
            
            # Check that the corresponding genome file exists
            genomeFile = os.path.join(inputDir, file.replace(".gff3", ".fasta"))
            if not os.path.exists(genomeFile):
                raise FileNotFoundError(f"Expected to find '{os.path.basename(genomeFile)}' in '{inputDir}' but could not")
            
            # Store the pairing
            filePairs.append([os.path.join(inputDir, file), genomeFile, suffixNum])
    
    # Filter down to only those that need to be generated
    filteredPrediction = []
    for gff3File, fastaFile, suffixNum in filePairs:
        # Derive output file names
        mrnaFileName = os.path.join(outputDir, f"{outputPrefix}{suffixNum}.mrna")
        cdsFileName = os.path.join(outputDir, f"{outputPrefix}{suffixNum}.cds")
        protFileName = os.path.join(outputDir, f"{outputPrefix}{suffixNum}.aa")
        sequenceFileNames = [mrnaFileName, cdsFileName, protFileName]
        
        # Generate files if they don't exist
        if not all([ os.path.exists(x) for x in sequenceFileNames ]) or \
        not all([ os.path.exists(x + ".ok") for x in sequenceFileNames ]):
            filteredPrediction.append([gff3File, fastaFile, [mrnaFileName, cdsFileName, protFileName]])
    
    # Extract annotations for files needing it
    if len(filteredPrediction) > 0:
        print(f"# Extracting sequences from GFF3/genome pair(s) to '{outputDir}' ...")
        
        for i in range(0, len(filteredPrediction), threads): # only process n (threads) files at a time
            processing = []
            for x in range(threads): # begin processing n files
                if i+x < len(filteredPrediction): # parent loop may excess if n > the number of files needing indexing
                    gff3File, fastaFile, (mrnaFileName, cdsFileName, protFileName) = filteredPrediction[i+x]
                    
                    extractionWorkerThread = GFF3ExtractionProcess(gff3File, fastaFile, mrnaFileName,
                                                                   cdsFileName, protFileName,
                                                                   isMicrobial, translationTable)
                    extractionWorkerThread.start()
                    processing.append(extractionWorkerThread)
            
            # Gather results
            for extractionWorkerThread in processing:
                extractionWorkerThread.join()
                extractionWorkerThread.check_errors()
