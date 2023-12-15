#! python3
# extract_other_sequences.py

# Script to take in the BINge working directory (which should have .nucl 
# files corresponding to mRNA/exon sequences in it) and generate .cds and
# .aa files for the GFF3s provided.

import os, argparse, sys
from Bio.Seq import Seq

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.fasta_handling import AnnotationExtractor

# Define functions
def validate_args(args):
    # Validate input directory location
    if not os.path.isdir(args.bingeDir):
        print(f'I am unable to locate the BINge working directory ({args.bingeDir})')
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    
    # Ensure that translationTable value is sensible
    if args.translationTable < 1:
        print('translationTable value must be greater than 1. Fix this and try again.')
        quit()
    elif args.translationTable > 31:
        print('translationTable value must be less than 31. Fix this and try again.')
        quit()

    # Validate output file location
    if os.path.isdir(args.outputDirectory):
        print(f'Output directory specified already exists ({args.outputDirectory})')
        print('This means I will try to resume a previously started run.')
        print('No files will be overwritten, so if something went wrong previously you ' + 
              'should fix that up before running this again.')
    elif not os.path.exists(args.outputDirectory):
        try:
            os.mkdir(args.outputDirectory)
            print(f'Output directory was created during argument validation')
        except Exception as e:
            print(f"An error occurred when trying to create the output directory '{args.outputDirectory}'")
            print("This probably means you've indicated a directory wherein the parent " + 
                  "directory does not already exist.")
            print("But, I'll show you the error below before this program exits.")
            print(e.message)
    else:
        print(f"Something other than a directory already exists at '{args.outputDirectory}'")
        print("Either move this, or specify a different -o value, then try again.")
        quit()

## Main
def main():
    # User input
    usage = """%(prog)s will produce CDS and protein translated sequences which should correspond
    to your annotation1.nucl, annotation2.nucl, ... etc., files generated as part of BINge's pipeline.
    These files may be helpful for use with unify_representatives.py or for any situation where you
    simply want those files.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="bingeDir",
                   required=True,
                   help="Input location where BINge was run")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output prefix for FASTA files")
    # Optional
    p.add_argument("--microbial", dest="isMicrobial",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if your BINge directory
                   points to GFF3 file(s) bacteria, archaea, or just any organism in
                   which the GFF3 does not contain mRNA and exon features; in this case,
                   I expect the GFF3 feature to have 'gene' and 'CDS' features.""",
                   default=False)
    p.add_argument("--translation", dest="translationTable",
                   required=False,
                   type=int, 
                   help="""Optionally specify the NCBI numeric genetic code to utilise for CDS
                   translation; this should be an integer from 1 to 31
                   (default == 1 i.e., Standard Code)""", default=1)
    args = p.parse_args()
    validate_args(args)
    
    # Locate subdirectory containing files
    gff3Dir = os.path.join(args.bingeDir, "gff3s")
    if not os.path.isdir(gff3Dir):
        print(f"I was expecting the directory '{gff3Dir}' to exist.")
        print("If your -i value is a BINge output directory, I should find 'gff3s' as a subdirectory.")
        print("Try to make this the case, or specify the correct -i location then try again.")
        quit()
    
    # Locate all GFF3/genome pairs
    filePairs = []
    for file in os.listdir(gff3Dir):
        if file.endswith(".gff3"):
            assert file.startswith("annotation"), \
                f"GFF3 file in '{gff3Dir}' has a different name than expected?"
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".gff3")[0]
            suffixNum = filePrefix.split("annotation")[1]
            assert suffixNum.isdigit(), \
                f"GFF3 file in '{gff3Dir}' does not have a number suffix?"
            
            # Check that the corresponding genome file exists
            genomeFile = f"genome{suffixNum}.fasta"
            assert os.path.exists(os.path.join(gff3Dir, genomeFile)), \
                f"Expected to find file '{genomeFile}' in directory '{gff3Dir}' but couldn't?"
            
            # Store the pairing
            filePairs.append([os.path.join(gff3Dir, file), os.path.join(gff3Dir, genomeFile), suffixNum])
    
    # Parse out CDS and protein sequences from each GFF3/genome pair
    warnOnce = False
    for gff3File, fastaFile, suffixNum in filePairs:
        cdsFileName = os.path.join(args.outputDirectory, f"annotations{suffixNum}.cds")
        aaFileName = os.path.join(args.outputDirectory, f"annotations{suffixNum}.aa")
        
        with open(cdsFileName, "w") as cdsOut, open(aaFileName, "w") as aaOut:
            seqGenerator = AnnotationExtractor(gff3File, fastaFile, args.isMicrobial)
            for mrnaID, exonSeq, cdsSeq in seqGenerator.iter_sequences():
                cdsOut.write(f">{mrnaID}\n{cdsSeq}\n")
                
                aaSeq = Seq(cdsSeq).translate(table=args.translationTable)
                aaOut.write(f">{mrnaID}\n{aaSeq}\n")
                
                if warnOnce == False:
                    if "*" in aaSeq[:-1]:
                        print(f"Warning: translation of '{mrnaID}' contains internal stop codons.")
                        print("Sequence may be flawed; further warnings will be suppressed.")
                        warnOnce = True
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
