import os, sys
from .validation import touch_ok
from .fasta_handling import AnnotationExtractor
from .thread_workers import BasicProcess

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_GFF3IO
from ZS_GFF3IO import GFF3

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
    def task(self, gff3FileIn, fastaFileIn, mrnaFileOut, cdsFileOut, protFileOut, isMicrobial):
        # Generate sequences
        with open(mrnaFileOut, "w") as mrnaOut, open(cdsFileOut, "w") as cdsOut, open(protFileOut, "w") as protOut:
            seqGenerator = AnnotationExtractor(gff3FileIn, fastaFileIn, isMicrobial)
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

def extract_annotations_from_gff3(locations, isMicrobial, threads):
    '''
    Will take the files within the 'gff3s' subdirectory of workingDirectory and
    produce sequence files from them.
    
    Parameters:
        locations -- a Locations object with attributes for directory locations.
        isMicrobial -- a boolean indicating whether the organism is a microbe and hence GFF3
                       has gene -> CDS features, rather than gene -> mRNA -> CDS/exon.
        threads -- an integer indicating the number of threads to use for processing.
    '''
    # Locate all GFF3/genome pairs
    filePairs = []
    for file in os.listdir(locations.gff3Dir):
        if file.endswith(".gff3"):
            if not file.startswith("annotation"):
                raise ValueError(f"'{file}' in '{locations.gff3Dir}' does not begin with 'annotation' as expected")
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".gff3")[0]
            suffixNum = filePrefix.split("annotation")[1]
            if not suffixNum.isdigit():
                raise ValueError(f"'{file}' in '{locations.gff3Dir}' does not have a number suffix as expected")
            
            # Check that the corresponding genome file exists
            genomeFile = os.path.join(locations.gff3Dir, f"genome{suffixNum}.fasta")
            if not os.path.exists(genomeFile):
                raise FileNotFoundError(f"Expected to find 'genome{suffixNum}.fasta' in '{locations.gff3Dir}' but could not")
            
            # Store the pairing
            filePairs.append([os.path.join(locations.gff3Dir, file), genomeFile, suffixNum])
    
    # Filter down to only those that need to be generated
    filteredPrediction = []
    for gff3File, fastaFile, suffixNum in filePairs:
        # Derive output file names
        mrnaFileName = os.path.join(locations.sequencesDir, f"annotations{suffixNum}.mrna")
        cdsFileName = os.path.join(locations.sequencesDir, f"annotations{suffixNum}.cds")
        protFileName = os.path.join(locations.sequencesDir, f"annotations{suffixNum}.aa")
        sequenceFileNames = [mrnaFileName, cdsFileName, protFileName]
        
        # Generate files if they don't exist
        if not all([ os.path.exists(x) for x in sequenceFileNames ]) or \
        not all([ os.path.exists(x + ".ok") for x in sequenceFileNames ]):
            filteredPrediction.append([gff3File, fastaFile, [mrnaFileName, cdsFileName, protFileName]])
    
    # Extract annotations for files needing it
    if len(filteredPrediction) > 0:
        print(f"# Generating sequences from GFF3/genome pair(s)...")
        
        for i in range(0, len(filteredPrediction), threads): # only process n (threads) files at a time
            processing = []
            for x in range(threads): # begin processing n files
                if i+x < len(filteredPrediction): # parent loop may excess if n > the number of files needing indexing
                    gff3File, fastaFile, (mrnaFileName, cdsFileName, protFileName) = filteredPrediction[i+x]
                    
                    extractionWorkerThread = GFF3ExtractionProcess(gff3File, fastaFile, mrnaFileName,
                                                                   cdsFileName, protFileName, isMicrobial)
                    extractionWorkerThread.start()
                    processing.append(extractionWorkerThread)
            
            # Gather results
            for extractionWorkerThread in processing:
                extractionWorkerThread.join()
                extractionWorkerThread.check_errors()
