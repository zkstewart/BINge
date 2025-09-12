import os, re, sys
from copy import deepcopy
from itertools import product
from Bio.Data import CodonTable, IUPACData

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from validation import touch_ok
from thread_workers import BasicProcess
from gff3 import GFF3Graph
from gff3fasta import FASTA

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
        # Parse genome sequence and annotation
        fasta = FASTA(fastaFileIn)
        gff3 = GFF3Graph(gff3FileIn)
        
        # Identify which feature types we will iterate through
        if isMicrobial:
            subfeature = ["gene"]
        else:
            subfeature = ["mRNA", "transcript"]
        
        # Generate sequences
        warnedOnce1 = False # trans-splicing warning
        warnedOnce2 = False # empty CDS warning
        warnedOnce3 = False # internal stop codon warning
        with open(mrnaFileOut, "w") as mrnaOut, open(cdsFileOut, "w") as cdsOut, open(protFileOut, "w") as protOut:
            for featureType in subfeature:
                if not featureType in gff3.ftypes:
                    continue
                
                for featureID in gff3.ftypes[featureType]:
                    feature = gff3[featureID]
                    
                    # Avoid trans-splicing
                    if feature.strand == "?":
                        if not warnedOnce1:
                            print(f"WARNING: Trans-splicing mixed strands (strand == '?') is not handled; '{featureID}' " + 
                                  f"and similar {featureType}'s will not be produced")
                            warnedOnce1 = True
                        continue
                    
                    # Avoid features without a CDS
                    if not hasattr(feature, "CDS"):
                        if not warnedOnce2:
                            print(f"WARNING: '{featureID}' lacks CDS annotation; it and similar {featureType}'s will not be part " + 
                                  "of the generated sequences")
                            warnedOnce2 = True
                        continue
                    
                    # Format the CDS and protein
                    cdsSequence = feature.as_gene_model(fasta, "CDS")
                    proteinSequence = cdsSequence.translate(translationTable)
                    
                    # Obtain the exon/mRNA depending on whether this is a microbial annotation
                    if isMicrobial:
                        exonSequence = cdsSequence
                    else:
                        exonSequence = feature.as_gene_model(fasta, "exon")
                    
                    # Write to file
                    mrnaOut.write(exonSequence.format())
                    cdsOut.write(cdsSequence.format())
                    protOut.write(proteinSequence.format())
        
        # Touch the .ok files to indicate completion
        touch_ok(mrnaFileOut)
        touch_ok(cdsFileOut)
        touch_ok(protFileOut)

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
