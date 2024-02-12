#! python3
# BINge.py
# BIN Genes for Expression analyses

# This script/program aims to provide an alternative to Corset
# for situations where a reference genome is available. However,
# the reference genome might be for a different subspecies than
# the one you're working with. Hence, it's not good enough to
# just map against the reference genome. You should map against a
# de novo transcriptome, but leverage the genomic information
# to group transcripts into genes. That's what this does.

import os, argparse, sys, pickle, platform, subprocess, json
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.bins import BinBundle
from modules.bin_handling import generate_bin_collections, populate_bin_collections
from modules.gmap_handling import setup_gmap_indices, auto_gmapping
from modules.clustering import cluster_unbinned_sequences
from modules.validation import validate_args, validate_fasta
from Various_scripts.Function_packages.ZS_Utility import convert_windows_to_wsl_path
from modules.fasta_handling import AnnotationExtractor, FastaCollection, \
    generate_sequence_length_index
HASHING_PARAMS = ["inputFiles", "genomeFiles",            # These hashes are the only ones 
                  "gmapIdentity", "clusterVoteThreshold"] # which behaviourally influence the 
                                                          # pre-external clustering
# Define functions
def symlinker(src, dst):
    '''
    Helper function to symlink files depending on host OS. For Linux, nothing special
    needs to occur. For Windows, developer mode is required. This function just acts as
    an exception handler.
    '''
    try:
        os.symlink(src, dst)
    except:
        if platform.system() == 'Windows':
            print("os.symlink is not working on your Windows computer.")
            print("This means developer mode is probably not activated.")
            print("Google 'enable windows developer mode' to see how to do this.")
            print("Until then, this program will exit since symlinks are required.")
        else:
            print(f"os.symlink failed for an unknown reason when linking '{src}' to '{dst}'.")
            print("This is unexpected and the program will exit.")
        quit()

def get_file_hash(originalFile):
    '''
    Creates a simple hash of the file name.
    
    Parameters:
        originalFile -- a string indicating the file name to hash.
    Returns:
        fileHash -- a string of the hash for the file name.
    '''
    return sha256(bytes(originalFile, 'utf-8')).hexdigest()

def write_file_hash(inputFile, linkedFile):
    '''
    Helper function to write a file with the hash of the input file. This is used
    to check if a file has changed across runs, and can help to unexpected results
    if that has occurred.
    
    Parameters:
        inputFile -- a string indicating the file name to hash.
        linkedFile -- a string indicating the file prefix to write the hash to.
    '''
    fileHash = get_file_hash(inputFile)
    open(f"{linkedFile}.{fileHash}", "w").close()

def check_file_hash(inputFile, linkedFile, errorIfNone=True):
    '''
    Helper function to check the hash of a file against the currently input file.
    If the hash is different, the program will exit.
    
    Parameters:
        inputFile -- a string indicating the file name to hash.
        linkedFile -- a string indicating the file prefix to write the hash to.
    '''
    SKIP_SUFFIXES = [ ".pkl", ".gmap", ".fai" ]
    
    fileHash = get_file_hash(inputFile)
    hashFileName = f"{linkedFile}.{fileHash}"
    
    # Find what hash files exist here
    hashFiles = [
        f
        for f in os.listdir(os.path.dirname(linkedFile))
        if f.startswith(os.path.basename(linkedFile))
        and not f == os.path.basename(linkedFile)
        and not any([ f.endswith(s) for s in SKIP_SUFFIXES ])
    ]
    hashSuffixes = [
        f.split(".")[-1]
        for f in hashFiles
    ]
    
    # Check if one, matching hash file exists
    if len(hashFiles) == 1 and os.path.exists(hashFileName):
        return # The file exists and the hash is the same; everything is okay
    
    # Check if one, non-matching file exists
    elif len(hashFiles) == 1:
        print(f"Expected to find file '{hashFileName}' but instead found '{hashFiles[0]}'...?")
        print(f"To fix this error, you should delete '{linkedFile}' and this hash file.")
        quit()
    
    # Check if more than one hash file exists
    elif len(hashFiles) > 1:
        print(f"Expected to find one hash file '{hashFileName}' but instead " + 
              f"found {len(hashFiles)} hash files...?")
        print(f"To fix this error, you should delete '{linkedFile}' and all of its hash files " +
              f"e.g., the files ending with values including {hashSuffixes}")
        quit()
    
    # Check if no hash files exist
    else:
        if errorIfNone:
            print(f"Expected to find file '{hashFileName}' but it doesn't exist...")
            print(f"To fix this error, you should delete '{linkedFile}' and rerun the program.")
            quit()
        else:
            write_file_hash(inputFile, linkedFile)

def check_file_exists(fileLocation):
    '''
    Helper function to check if a symbolic link/file exists. Using os.path.exists() doesn't
    show an existing symlink on Windows if it was made through WSL. So we need to
    think outside the box to check if it exists or not.
    '''
    if platform.system() != 'Windows':
        return os.path.exists(fileLocation)
    else:
        fileLocation = convert_windows_to_wsl_path(fileLocation)
        cmd = ["wsl", "~", "-e", "ls", fileLocation]
        run_file_exists = subprocess.Popen(cmd, shell = True,
                                           stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        existsout, existserr = run_file_exists.communicate()
        if "No such file or directory" in existserr.decode("utf-8"):
            return False
        elif existsout.decode("utf-8").rstrip("\r\n ") == fileLocation:
            return True
        else:
            raise Exception(("symlink_exists encountered an unhandled situation; have a look " +
                             "at the stdout '" + existsout.decode("utf-8") + "' and stderr '" + 
                             existserr.decode("utf-8") + "' to make sense of this."))

def _setup_error_helper(fileDir, filesInDir, numFiles, fileType):
    '''
    Simple function to be called by setup_working_directory to raise an informative
    error message if we detect that an analysis is being re-run inappropriately.
    '''
    if len(filesInDir) > numFiles:
        print(f"Expected to find {numFiles} {fileType} files at '{fileDir}' but instead found {len(filesInDir)}...")
        print("Are you trying to run an analysis in an existing dir but removing files from your input arguments?")
        print("This isn't supported; please use a fresh directory for each analysis.")
        quit()

def setup_working_directory(fileNames, genomeFiles, workingDirectory):
    '''
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        fileNames -- a list of strings pointing to FASTA or GFF3 files.
        genomeFiles -- a list of strings pointing to genome FASTAs to align against.
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
    '''
    # Create subdirectory for files (if not already existing)
    gff3Dir = os.path.join(workingDirectory, "gff3s")
    genomesDir = os.path.join(workingDirectory, "genomes")
    os.makedirs(gff3Dir, exist_ok=True)
    os.makedirs(genomesDir, exist_ok=True)
    
    # Link to the -i fileNames values
    numGFF3s = 0
    numFASTAs = 0
    for file in fileNames:
        # Handle GFF3:FASTA pairs
        if "," in file:
            gff3, fasta = file.split(",")
            numGFF3s += 1
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                print(f"-i value '{file}' value after the ',' is not a FASTA file")
                print("Make sure you specify the file order as GFF3:FASTA then try again.")
                quit()
            
            # Symlink files to GFF3s subdirectory if not aleady existing
            linkedGFF3 = os.path.join(gff3Dir, f"annotation{numGFF3s}.gff3")
            linkedFASTA = os.path.join(gff3Dir, f"genome{numGFF3s}.fasta")
            
            if not check_file_exists(linkedGFF3):
                symlinker(gff3, linkedGFF3)
                write_file_hash(gff3, linkedGFF3)
            else:
                check_file_hash(gff3, linkedGFF3)
            
            if not check_file_exists(linkedFASTA):
                symlinker(fasta, linkedFASTA)
                write_file_hash(fasta, linkedFASTA)
            else:
                check_file_hash(fasta, linkedFASTA)
        
        # Handle plain FASTA files
        else:
            numFASTAs += 1
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(file)
            if not isFASTA:
                print(f"-i value '{file}' is not a FASTA file")
                print("Make sure you specify the right file and/or location then try again.")
                quit()
            
            # Symlink to main working directory if not already existing
            linkedTranscriptome = os.path.join(workingDirectory, f"transcriptome{numFASTAs}.nucl")
            
            if not check_file_exists(linkedTranscriptome):
                symlinker(file, linkedTranscriptome)
                write_file_hash(file, linkedTranscriptome)
            else:
                check_file_hash(file, linkedTranscriptome)
    
    # Check that GFF3 and FASTAs do not have excess
    gff3sInDir = [ f for f in os.listdir(gff3Dir) if f.endswith(".gff3") ]
    fastasInDir = [ f for f in os.listdir(gff3Dir) if f.endswith(".fasta") ]
    _setup_error_helper(gff3Dir, gff3sInDir, numGFF3s, "GFF3")
    _setup_error_helper(gff3Dir, fastasInDir, numFASTAs, "FASTA")
    
    # Link to the -g genomeFiles values
    numGFF3s = 0
    numGenomes = 0
    for file in genomeFiles:
        numGenomes += 1
        # Handle GFF3:FASTA pairs
        if "," in file:
            numGFF3s += 1
            gff3, fasta = file.split(",")
            
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(fasta)
            if not isFASTA:
                print(f"-g value '{fasta}' value after the ',' is not a FASTA file")
                print("Make sure you specify the file order as GFF3:FASTA then try again.")
                quit()
            
            # Symlink files to genomes subdirectory if not aleady existing
            linkedGFF3 = os.path.join(genomesDir, f"annotation{numGenomes}.gff3")
            linkedFASTA = os.path.join(genomesDir, f"genome{numGenomes}.fasta")
            
            if not check_file_exists(linkedGFF3):
                symlinker(gff3, linkedGFF3)
                write_file_hash(gff3, linkedGFF3)
            else:
                check_file_hash(gff3, linkedGFF3)
            
            if not check_file_exists(linkedFASTA):
                symlinker(fasta, linkedFASTA)
                write_file_hash(fasta, linkedFASTA)
            else:
                check_file_hash(fasta, linkedFASTA)
            
            # Index the genome's contig lengths if not already done
            generate_sequence_length_index(linkedFASTA)
        
        # Handle plain genome files
        else:
            # Check that FASTA is a FASTA
            isFASTA = validate_fasta(file)
            if not isFASTA:
                print(f"-g value '{file}' is not a FASTA file")
                print("Make sure you specify the right file and/or location then try again.")
                quit()
            
            # Symlink to main working directory if not already existing
            linkedFASTA = os.path.join(genomesDir, f"genome{numGenomes}.fasta")
            
            if not check_file_exists(linkedFASTA):
                symlinker(file, linkedFASTA)
                write_file_hash(file, linkedFASTA)
            else:
                check_file_hash(file, linkedFASTA)
            
            # Index the genome's contig lengths if not already done
            generate_sequence_length_index(linkedFASTA)
    
    # Check that genomes do not have excess
    gff3sInDir = [ f for f in os.listdir(genomesDir) if f.endswith(".gff3") ]
    genomesInDir = [ f for f in os.listdir(genomesDir) if f.endswith(".fasta") ]
    _setup_error_helper(genomesDir, gff3sInDir, numGFF3s, "GFF3")
    _setup_error_helper(genomesDir, genomesInDir, numGenomes, "FASTA")

def load_param_cache(workingDirectory):
    '''
    Loads the parameter cache file from the working directory, if it exists,
    as a dictionary with JSON parsing.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
    '''
    # Parse any existing param cache file
    paramCacheFile = os.path.join(workingDirectory, "param_cache.json")
    if os.path.exists(paramCacheFile):
        try:
            with open(paramCacheFile, "r") as fileIn:
                paramsDict = json.load(fileIn)
        except:
            raise Exception((f"'{paramCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it so I can reinitialise one."))
    else:
        paramsDict = {}
    
    return paramsDict

def setup_param_cache(args, paramHash):
    '''
    Writes the parameter hash to the param cache file in the working directory, appending
    the parameters to the cache if one already exists.
    
    Parameters:
        args -- the argparse object of BINge called through the main function.
        paramHash -- a string of the hash for these parameters to store in the param cache.
    '''
    # Parse any existing param cache file
    paramsDict = load_param_cache(args.outputDirectory)
    
    # Add this program run to the paramsDict cache (if needed)
    paramsDict[paramHash] = {
        param : args.__dict__[param]
        for param in HASHING_PARAMS
    }
    
    # Write updated param cache to file
    with open(os.path.join(args.outputDirectory, "param_cache.json"), "w") as fileOut:
        json.dump(paramsDict, fileOut)

def setup_sequences(workingDirectory, isMicrobial=False):
    '''
    Will take the files within the 'gff3s' subdirectory of workingDirectory and
    produce sequence files from them.
    
    Given a mix of FASTA and/or GFF3 files, this function will symlink FASTA files
    and generate FASTAs from the GFF3 files at the indicated working directory location.
    
    Parameters:
        workingDirectory -- a string indicating an existing directory to symlink and/or
                            write FASTAs to.
        isMicrobial -- a boolean indicating whether the organism is a microbe and hence GFF3
                       has gene -> CDS features, rather than gene -> mRNA -> CDS/exon.
    '''
    # Locate subdirectory containing files
    gff3Dir = os.path.join(workingDirectory, "gff3s")
    assert os.path.isdir(gff3Dir), \
        f"setup_sequences failed because '{gff3Dir}' doesn't exist somehow?"
    
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
            genomeFile = os.path.join(gff3Dir, f"genome{suffixNum}.fasta")
            assert check_file_exists(genomeFile), \
                f"Expected to find file 'genome{suffixNum}.fasta' at '{gff3Dir}' but couldn't?"
            
            # Store the pairing
            filePairs.append([os.path.join(gff3Dir, file), genomeFile, suffixNum])
    
    # Parse out mRNA sequences from each GFF3/genome pair
    for gff3File, fastaFile, suffixNum in filePairs:
        sequenceFileName = os.path.join(workingDirectory, f"annotations{suffixNum}.nucl")
        
        # Generate file if it doesn't exist
        if not check_file_exists(sequenceFileName):
            try:
                with open(sequenceFileName, "w") as fileOut:
                    seqGenerator = AnnotationExtractor(gff3File, fastaFile, isMicrobial)
                    for mrnaID, exonSeq, cdsSeq in seqGenerator.iter_sequences():
                        fileOut.write(f">{mrnaID}\n{cdsSeq}\n")
            except Exception as e:
                if check_file_exists(sequenceFileName):
                    os.unlink(sequenceFileName)
                raise e

def find_missing_sequence_id(binCollectionList, transcriptRecords):
    '''
    Compares one or more BinCollection objects against the transcript sequences
    to see if any sequences indicated in the BinCollection do not exist in
    transcriptRecords. If any are found to be missing, the first one will be returned
    as an exemplar.
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        missingSeqID -- a string of the first sequence encountered in our BinCollections
                        that could not be found in transcriptRecords.
    '''
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            for seqID in bin.ids:
                if seqID not in transcriptRecords:
                    return seqID
    return None

def get_unbinned_sequence_ids(clusterDict, eliminatedIDs, transcriptRecords):
    '''
    Compares one or more BinBundle objects against the transcript sequences
    to see if any sequences indicated in transcriptRecords do not exist in
    any Bins.
    
    Parameters:
        clusterDict -- a dictionary with structure like:
                       {
                           0 : [ "seq1", "seq2", "seq3" ],
                           1 : [ "seq4", "seq5", "seq6" ],
                           ...
                        }
        eliminatedIDs -- a set containing strings of sequence IDs which are not to be
                         considered for clustering.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        unbinnedIDs -- a set containing string value for sequence IDs that were not
                       binned by BINge's main clustering process.
    '''
    binnedIDs = []
    for seqIDs in clusterDict.values():
        binnedIDs.extend(seqIDs)
    binnedIDs = set(binnedIDs)
    binnedIDs = binnedIDs.union(eliminatedIDs)
    
    unbinnedIDs = set()
    for record in transcriptRecords:
        if record.name not in binnedIDs:
            unbinnedIDs.add(record.name)
    
    return unbinnedIDs

def get_parameters_hash(args):
    '''
    Function to receive the arguments associated with BINge.py and generate a hash
    unique to this run. This hash can be used for pickle persistance of data which
    will enable the program to resume from the point that BINge's operations finish
    and the external clustering begins.
    
    Parameters:
        args -- the Argparse object associated with the main BINge.py function.
    Returns:
        paramHash -- a sha256 hash string of the parameters
    '''
    strForHash = ""
    for param in HASHING_PARAMS:
        strForHash += str(args.__dict__[param])
    paramHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return paramHash

## Main
def main():
    showHiddenArgs = '--help-long' in sys.argv
    
    # User input
    usageShort = """Quick notes for the use of %(prog)s:
    1) if providing GFF3 as input, only protein-coding mRNA features will be handled; input
    FASTA file(s) are similarly expected to be nucleotide CDS or mRNA sequences.
    2) when providing GFF3 as input, you must also specify the genome FASTA it's associated
    with; you can do that with the format 'file.gff3,file.fasta'
    3) provide as many genome files as are relevant to your analysis, alongside their
    annotations (if available and high quality) in the same style e.g., 
    'annotation.gff3,genome.fasta'
    4) Unbinned sequences will be cascade clustered using MMseqs2 by default which is
    recommended; you can choose Linclust or CD-HIT if desired.
    5) Many parameters are hidden in this short help format since their defaults are adequate;
    specify --help-long to see information for those options.
    """
    
    usageLong = """%(prog)s (BIN Genes for Expression analyses) is a program which clusters
    (or "bins") de novo-assembled transcripts together with reference genome models if they
    are available. This might be necessary when working with multiple subspecies that are
    expected to diverge only slightly from the reference organism. By binning like this, each
    subspecies can have its gene counts compared fairly during DGE, and some of the pitfalls
    of other approaches e.g., CD-HIT are avoided.
    ###
    This program requires you to have GMAP and either MMseqs2 or CD-HIT available on your
    system. It will automatically make use of these softwares as part of its operations.
    ###
    In short, BINge operates by aligning all input files (-i) against all genomes (-g). In
    doing so, isoforms will cluster by their common best alignment positions, and synteny
    will be implicitly inferred based on where each species' best alignment positions are across
    the input genomes. This biologically minded clustering approach is likely to outperform
    methods which solely look at sequence identity.
    ###
    Note 1: You should try to provide a genome file for each species you've used in the -i
    argument or, minimally, provide the most closely related genome for each one. This isn't
    strictly necessary but you should try to do this wherever possible.
    ###
    Note 2: If you provide a GFF3 in the -i argument, you must indicate the genome FASTA it
    is associated with by linking the files together as a single string. For example you might
    provide inputs like '-i transcriptome.fasta annotation.gff3,reference.fasta'; as indicated,
    a standlone FASTA is okay, but a GFF3 must be paired with its sequence file via the 
    ',' character and in the order indicated i.e., annotation GFF3 then reference FASTA.
    ###
    Note 3: Similarly, for genome FASTAs given in the -g argument, you should either provide the
    genome on its own, or alongside its annotation GFF3. For example here you might provide inputs
    like '-g genome1.fasta annotation2.gff3,genome2.fasta'; by providing a GFF3 alongside the
    genome it will pre-seed the bins along this genome based on the annotation. If the annotation
    is of a reasonable standard this is expected to make BINge perform better.
    ###
    Note 4: Sequences which do not align against the genome are considered to be "unbinned".
    These will be clustered with MMseqs2 cascaded clustering by default, which is the recommended
    choice. You can use Linclust (also okay, trades some accuracy for some speed) or CD-HIT
    (potentially very slow and possibly least accurate) if wanted.
    ###
    Note 5: The --gmapIdentity parameter should be set in the range of 0.90 to 0.99 depending
    on the evolutionary distance between your input files and the genomes you're aligning
    against. For same species, use 0.98 or 0.99. If you're aligning against a different species
    in the same genus, use 0.95. If you're aligning against a different genus, consider 0.90.
    ###
    Note 6: The --clusterVoteThreshold, based on objective evidence, should be set to 0.5 or
    0.66. Using a lower value will give a tighter clustering (fewer clusters) and a higher
    value will give looser clustering (more clusters). A value of 0.66 is likely to be more
    biologically correct based on objective evaluation, but the difference is very marginal
    in terms of biological correctness and number of clusters (e.g., you may find up to 3%
    more clusters in 0.66 relative to 0.5). Generally, just stick to the default and you'll
    be alright.
    ###
    Note 7: You're seeing the --help-long format of this message, which means you may want to
    configure the way clustering of unbinned sequences works. Behavioural parameters of several
    features can be tuned here, but the defaults are expected to work most of the time.
    The main exception is CD-HIT's memory utilisation, which probably should be set depending
    on what you have available if you choose to use CD-HIT. tldr; change these only if you know
    what you're doing since they have been set to defaults which are likely to be optimal for
    most use cases.
    """
    
    p = argparse.ArgumentParser(description=usageLong if showHiddenArgs else usageShort)
    # Required
    p.add_argument("-i", dest="inputFiles",
                   required=True,
                   nargs="+",
                   help="""Input transcriptome FASTA(s) and/or annotation GFF3(s) paired to
                   their genome file with ',' separator""")
    p.add_argument("-g", dest="genomeFiles",
                   required=True,
                   nargs="+",
                   help="Input genome FASTA(s) to align against")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory for intermediate and final results")
    # Optional - BINge
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many threads to run when multithreading
                   is available (default==1)""",
                   default=1)
    p.add_argument("--clusterVoteThreshold", dest="clusterVoteThreshold",
                   required=False,
                   type=float,
                   help="""Optionally, specify the clustering vote threshold used when
                   clustering bins based on sequence ID co-occurrence (default == 0.66).
                   A higher value gives more clusters and a lower value gives fewer
                   clusters; if you want to be more biologically correct, use 0.66, but
                   if you want fewer gene clusters, use 0.5.
                   """
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.66)
    # Optional - GMAP
    p.add_argument("--gmapDir", dest="gmapDir",
                   required=False,
                   help="""If GMAP is not discoverable in your PATH, specify the directory
                   containing the mmseqs executable""")
    p.add_argument("--gmapIdentity", dest="gmapIdentity",
                   required=False,
                   type=float,
                   help="""Optionally, specify the identity threshold for accepting a GMAP
                   alignment (default==0.95); note that this value operates independently
                   of --identity and its strictness should depend on the largest evolutionary
                   distance you have between a file given to -i and a genome given to -g e.g.,
                   this should be strict for same species only alignment, less strict for
                   same genus alignment, and least strict for different genus alignments"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.95)
    # Optional - program behavioural controls
    p.add_argument("--microbial", dest="isMicrobial",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you are providing a GFF3
                   file from a bacteria, archaea, or just any organism in which the GFF3
                   does not contain mRNA and exon features; in this case, I expect the GFF3
                   feature to have 'gene' and 'CDS' features."""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=False)
    p.add_argument("--clusterer", dest="unbinnedClusterer",
                   required=False,
                   choices=["mmseqs-cascade", "mmseqs-linclust", "cd-hit"],
                   help="""Specify which algorithm to use for clustering of unbinned sequences
                   (default=='mmseqs-cascade')""",
                   default="mmseqs-cascade")
    p.add_argument("--mmseqs", dest="mmseqsDir",
                   required=False,
                   help="""If using MMseqs2-based clustering, specify the directory containing
                   the mmseqs executable""")
    p.add_argument("--cdhit", dest="cdhitDir",
                   required=False,
                   help="""If using CD-HIT clustering, specify the directory containing
                   the cd-hit-est executable""")
    p.add_argument("--identity", dest="identity",
                   required=False,
                   type=float,
                   help="""ALL CLUSTERERS: Specify the identity threshold for clustering
                   (default==0.98); this value should be strict unless you are clustering
                   multiple species' together for a DGE analysis"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.98)
    # Optional - MMseqs2
    p.add_argument("--tmpDir", dest="tmpDir",
                   required=False,
                   help="""MMSEQS: Specify the tmpDir for MMseqs2 running; default='mms2_tmp'
                   in your current working directory"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="mms2_tmp")
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="MMSEQS: Specify the evalue threshold for clustering (default==1e-3)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=1e-3)
    p.add_argument("--coverage", dest="coverage",
                   required=False,
                   type=float,
                   help="MMSEQS: Specify the coverage ratio for clustering (default==0.4)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--mode", dest="mode",
                   required=False,
                   choices=["set-cover", "connected-component", "greedy"],
                   help="MMSEQS: Specify the clustering mode (default=='connected-component')"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="connected-component")
    p.add_argument("--sensitivity", dest="sensitivity",
                   required=False,
                   choices=["4","5","5.7","6","7","7.5"],
                   help="MMSEQS-CASCADE: Specify the sensitivity value (default==5.7)"
                   if showHiddenArgs else argparse.SUPPRESS,
                   default="5.7")
    p.add_argument("--steps", dest="steps",
                   required=False,
                   type=int,
                   help="""MMSEQS-CASCADE: Specify the number of cascaded clustering steps 
                   (default==3)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=3)
    # Optional - CD-HIT
    p.add_argument("--cdhit_shortcov", dest="cdhitShortCov",
                   required=False,
                   type=float,
                   help="""CDHIT: Specify what -aS parameter to provide
                   CD-HIT (default==0.4)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--cdhit_longcov", dest="cdhitLongCov",
                   required=False,
                   type=float,
                   help="""CDHIT: Specify what -aL parameter to provide
                   CD-HIT (default==0.4)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=0.4)
    p.add_argument("--cdhit_mem", dest="cdhitMem",
                   required=False,
                   type=int,
                   help="""CDHIT: Specify how many megabytes of memory to
                   provide CD-HIT (default==6000)"""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=6000)
    # Help controller and meta arguments
    p.add_argument("--debug", dest="debug",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you want to generate detailed
                   logging information along the way to help with debugging."""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=False)
    p.add_argument("--help-long", dest="help-long",
                   action="help",
                   help="""Show all options, including those that are not
                   recommended to be changed"""
                   if not showHiddenArgs else argparse.SUPPRESS)
    
    args = p.parse_args()
    validate_args(args)
    paramHash = get_parameters_hash(args)
    print(f"# Hash associated with this run is {paramHash}")
    
    # Setup sequence working directory for analysis
    setup_working_directory(args.inputFiles, args.genomeFiles, args.outputDirectory)
    setup_param_cache(args, paramHash)
    
    # Figure out if we've already run BINge here before and exit if so
    outputFileName = os.path.join(args.outputDirectory, f"BINge_clustering_result.{paramHash}.tsv")
    if os.path.exists(outputFileName):
        print(f"A BINge output file already exists at '{outputFileName}'")
        print("This program will not overwrite an existing file.")
        print("If you want to resume an existing run, make sure to move/rename/delete this file first.")
        print("Fix this before trying again.")
        quit()
    
    # Extract mRNAs from any input GFF3 annotations
    setup_sequences(args.outputDirectory, args.isMicrobial)
    
    # Establish GMAP indexes
    setup_gmap_indices(args.outputDirectory, args.gmapDir, args.threads)
    
    # Perform GMAP mapping
    gmapFiles = auto_gmapping(args.outputDirectory, args.gmapDir, args.threads)
    
    # Figure out what our pickle file is called
    pickleFile = os.path.join(args.outputDirectory, f"{paramHash}.binge.pkl")
    
    # Either load a pickle generated by previous BINge run ...
    if os.path.isfile(pickleFile) or os.path.islink(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            binBundle = pickle.load(pickleIn)
        print("Note: an existing pickle file with a matching parameters hash is present " + 
              f"i.e., '{pickleFile}'.")
        print("I will load that in now and resume program operation.")
    
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        print(f"{pickleFile} already exists, but is not a file?")
        print("BINge expects this to be a file which it can read or write to.")
        print("Something weird is happening, so I will exit the program now.")
        print(f"Move whatever is at the location of '{pickleFile}' then try again.")
        quit()
    
    # ... or begin pre-external clustering BINge
    else:
        # Set up a bin collection structure for each genome
        collectionList = generate_bin_collections(args.outputDirectory, args.threads, args.isMicrobial)
        if args.debug:
            print(f"# Generated a list with {len(collectionList)} collections")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} contains {len(_cl)} bins")
        
        # Parse GMAP alignments into our bin collection with multiple threads
        collectionList = populate_bin_collections(args.outputDirectory,
                                                  collectionList, gmapFiles,
                                                  args.threads, args.gmapIdentity)
        if args.debug:
            print(f"# Populated collections based on GMAP alignments")
            for index, _cl in enumerate(collectionList):
                print(f"# Collection #{index+1} now contains {len(_cl)} bins")
        
        # Convert collections into a bundle
        binBundle = BinBundle.create_from_multiple_collections(collectionList)
        
        # Cluster bundles across and within genomes
        clusterDict, eliminations = binBundle.cluster_by_cooccurrence(args.clusterVoteThreshold)
        if args.debug:
            print(f"# Clustered bundles (across and within genomes) based on ID occurrence")
            print(f"# Cluster dictionary contains {len(clusterDict)} clusters")
        
        # Write pickle file for potential resuming of program
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(clusterDict, pickleOut)
    
    # Write binned clusters to file
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines
        for clusterNum, seqIDs in clusterDict.items():
            for seqID in seqIDs:
                fileOut.write(f"{clusterNum}\t{seqID}\tbinned\n")
    
    # Cluster remaining unbinned sequences
    transcriptRecords = FastaCollection([
        os.path.join(args.outputDirectory, f)
        for f in os.listdir(args.outputDirectory)
        if f.endswith(".nucl")
    ])
    unbinnedIDs = get_unbinned_sequence_ids(clusterDict, eliminations, transcriptRecords)
    if args.debug:
        print(f"# There are {len(unbinnedIDs)} unbinned sequences for external clustering")
    
    if len(unbinnedIDs) == 1:
        unbinnedClusterDict = { 0: list(unbinnedIDs)[0] }
    elif len(unbinnedIDs) > 0:
        unbinnedClusterDict = cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args)
    else:
        unbinnedClusterDict = {} # blank to append nothing to output file
    
    # Write output of clustering to file
    numClusters = len(clusterDict)
    with open(outputFileName, "a") as fileOut:
        for clusterNum, clusterIDs in unbinnedClusterDict.items():
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum+numClusters+1}\t{seqID}\tunbinned\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
