import os, re, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from validation import touch_ok
from thread_workers import BasicProcess
from gff3 import GFF3Graph
from gff3fasta import FASTA

def gff3_to_fasta(gff3FileIn, fastaFileIn, mrnaFileOut, cdsFileOut, protFileOut, isMicrobial, translationTable):
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
                        print(f"WARNING: '{featureID}' in '{gff3FileIn}' has trans-splicing of mixed strands (strand == '?') " + 
                              f"which is not handled by BINge; it and similar {featureType}'s will not be produced")
                        warnedOnce1 = True
                    continue
                
                # Avoid features without a CDS
                if not hasattr(feature, "CDS"):
                    if not warnedOnce2:
                        print(f"WARNING: '{featureID}' in '{gff3FileIn}' lacks CDS annotation; it and similar " + 
                              f"{featureType}'s will not be produced")
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
