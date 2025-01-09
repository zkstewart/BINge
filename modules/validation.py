import os, distutils.spawn, platform, subprocess, sys
from pathlib import Path

from .locations import Locations

# Argument validations
def _validate_gmap(args):
    if args.gmapDir == None:
        gmap = distutils.spawn.find_executable("gmap")
        gmap_build = distutils.spawn.find_executable("gmap_build")
        
        # Raise error if we failed to find the executables
        if gmap == None or gmap_build == None:
            raise FileNotFoundError("--gmapDir wasn't specified, and either 'gmap' or 'gmap_build' are " + 
                                    "missing from your PATH")
        else:
            gmapDir = os.path.dirname(gmap)
            gmap_buildDir = os.path.dirname(gmap_build)
            
            # Raise errors if they're not in the same spot
            if gmapDir != gmap_buildDir:
                raise ValueError("--gmapDir wasn't specified, and I found 'gmap' or 'gmap_build' in " + 
                                 "different locations in your PATH. Make sure they're in the same " +
                                 "directory and try again.")
            
            # Set the value so we can use it later
            args.gmapDir = gmapDir
    else:
        if not os.path.isdir(args.gmapDir):
            raise FileNotFoundError(f"Unable to locate the GMAP directory '{args.gmapDir}'")

def initialise_working_directory(workingDirectory):
    # Validate working directory
    workingDirectory = os.path.abspath(workingDirectory)
    if os.path.isdir(workingDirectory):
        print("Working directory already exists; will attempt to resume a previous initialisation...")
    elif not os.path.exists(workingDirectory):
        try:
            os.mkdir(workingDirectory)
            print(f'Working directory was created during argument validation')
        except Exception as e:
            print(f"An error occurred when trying to create the working directory '{workingDirectory}'")
            print("This probably means you've indicated a directory wherein the parent " + 
                  "directory does not already exist.")
            print("I'll show you the error below before this program exits.")
            raise Exception(e.message)
    else:
        raise ValueError(f"Something other than a directory already exists at '{workingDirectory}'. " +
                         "Please move this, or specify a different -d value, then try again.")

def validate_init_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory) # performs validation implicitly
    
    # Validate -ig file locations
    for inputArgument in args.inputGff3Files:
        if not inputArgument.count(",") == 1:
            raise ValueError(f"-ig value '{inputArgument}' must have 1 comma in it separating GFF3,genome files")
        
        for inputFile in inputArgument.split(","):
            if not os.path.isfile(inputFile):
                raise ValueError(f"Unable to locate the file '{inputFile})' from the -ig argument '{inputArgument}'")
    
    # Validate -ix file locations
    for inputArgument in args.inputTxomeFiles:
        if not inputArgument.count(",") in [0, 2]:
            raise ValueError(f"-ix value '{inputArgument}' must have 0 or two commas in it (commas would separate " +
                             "mRNA,CDS,protein files)")
        
        for inputFile in inputArgument.split(","):
            if not os.path.isfile(inputFile):
                raise ValueError(f"Unable to locate the file '{inputFile})' from the -ix argument '{inputArgument}'")
    
    # Validate -t file locations
    for inputArgument in args.targetGenomeFiles:
        if not inputArgument.count(",") <= 1:
            raise ValueError(f"-t value '{inputArgument}' must have 0 or one comma in it (comma would separate GFF3,genome files)")
        
        for inputFile in inputArgument.split(","):
            if not os.path.isfile(inputFile):
                raise FileNotFoundError(f"Unable to locate the -t input file '{inputFile}'")
    
    # Validate numeric parameters
    if args.threads < 1:
        raise ValueError("--threads should be given a value >= 1")
    
    # Validate GMAP location
    _validate_gmap(args)
    
    return locations

def validate_cluster_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Validate numeric BINge parameters
    if args.threads < 1:
        raise ValueError("--threads should be given a value >= 1")
    if not 0.0 < args.identity <= 1.0:
        raise ValueError("--identity should be given a value greater than zero, and equal to " + 
                         "or less than 1")
    if not 0.0 < args.gmapIdentity <= 1.0:
        raise ValueError("--gmapIdentity should be given a value greater than zero, and equal to " + 
                         "or less than 1")
    if not 0.0 < args.clusterVoteThreshold <= 1.0:
        raise ValueError("--clusterVoteThreshold should be given a value greater than zero, and equal to " + 
                         "or less than 1")
    
    # Validate GMAP location
    _validate_gmap(args)
    
    # Validate optional MMseqs2 parameters
    if args.unbinnedClusterer in ["mmseqs-cascade", "mmseqs-linclust"]:
        # Validate MMseqs2 location (if applicable)
        if args.mmseqsDir == None:
            mmseqs = distutils.spawn.find_executable("mmseqs")
            if mmseqs == None:
                raise FileNotFoundError("--mmseqsDir wasn't specified, and 'mmseqs' is missing from your PATH")
            args.mmseqsDir = os.path.dirname(mmseqs)
        
        # Validate remaining parameters
        if not os.path.isdir(args.mmseqsDir):
            raise FileNotFoundError(f"Unable to locate the MMseqs2 directory '{args.mmseqsDir}'")
        if args.mmseqsEvalue < 0:
            raise ValueError("--mmseqs_evalue must be greater than or equal to 0")
        if not 0 <= args.mmseqsCoverage <= 1.0:
            raise ValueError("--mmseqs_cov must be a float in the range 0.0 -> 1.0")
        "--mode is controlled by argparse choices"
        
        # Validate "MMS-CASCADE" parameters
        if args.mmseqsSensitivity in ["5.7", "7.5"]:
            args.mmseqsSensitivity = float(args.mmseqsSensitivity)
        else:
            args.mmseqsSensitivity = int(args.mmseqsSensitivity)
        
        if args.mmseqsSteps < 1:
            raise ValueError("--mmseqs_steps must be greater than or equal to 1")
    
    # Validate optional CD-HIT parameters
    if args.unbinnedClusterer == "cd-hit":
        # Validate CD-HIT location (if applicable)
        if args.cdhitDir == None:
            cdhitest = distutils.spawn.find_executable("cd-hit-est")
            if cdhitest == None:
                raise FileNotFoundError("--cdhitDir wasn't specified, and 'cd-hit-est' is missing from your PATH")
            args.cdhitDir = os.path.dirname(cdhitest)
        
        # Validate remaining parameters
        if not os.path.isdir(args.cdhitDir):
            raise FileNotFoundError(f"Unable to locate the CD-HIT directory '{args.cdhitDir}'")
        if not 0.0 <= args.cdhitShortCov <= 1.0:
            raise ValueError("--cdhit_shortcov should be given a value in the range of 0 -> 1 (inclusive)")
        if not 0.0 <= args.cdhitLongCov <= 1.0:
            raise ValueError("--cdhit_longcov should be given a value in the range of 0 -> 1 (inclusive)")
        if not args.cdhitMem >= 100:
            raise ValueError("--cdhit_mem should be given a value at least greater than 100 (megabytes)")
    
    return locations

def validate_view_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Validate that viewed directory exists
    locations.runName = args.analysisFolder
    args.runDir = locations.resolve_runName(locations.analysisDir) # performs validation implicitly
    
    return locations

def validate_blast_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Validate that sequences directory exists
    sequencesDir = os.path.join(args.workingDirectory, "sequences")
    if not os.path.isdir(sequencesDir):
        raise FileNotFoundError(f"Unable to locate 'sequences' within '{sequencesDir}'; " +
                                "have you run the initialisation step yet?")
    
    # Locate and validate that sequences exist in the sequences directory
    sequenceSuffix = ".cds" if args.sequenceType == "nucleotide" else ".aa"
    args.sequenceFiles = validate_sequence_files(sequencesDir, sequenceSuffix)
    
    # Validate that the BLAST database exists
    args.targetFile = os.path.abspath(args.targetFile)
    if not os.path.isfile(args.targetFile):
        raise FileNotFoundError(f"Unable to locate the -t FASTA file '{args.targetFile}'")
    
    # Validate that MMseqs2 executable is locateable
    if args.mms2Exe == None:
        mms2 = distutils.spawn.find_executable("mmseqs")
        if mms2 == None:
            raise FileNotFoundError("--mms2 wasn't specified, and 'mmseqs' is missing from your PATH")
        args.mms2Exe = mms2
    else:
        if not os.path.isfile(args.mms2Exe):
            raise FileNotFoundError(f"Unable to locate the MMseqs2 executable '{args.mms2Exe}'")
    args.mms2Exe = os.path.abspath(args.mms2Exe)
    
    # Validate numeric parameters
    if args.threads < 1:
        raise ValueError("--threads should be given a value >= 1")
    
    return locations

def validate_salmon_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Validate that sequences directory exists
    sequencesDir = os.path.join(args.workingDirectory, "sequences")
    if not os.path.isdir(sequencesDir):
        raise FileNotFoundError(f"Unable to locate 'sequences' within '{sequencesDir}'; " +
                                "have you run the initialisation step yet?")
    
    # Validate that reads directory exists
    for rDir in args.readsDir:
        if not os.path.isdir(rDir):
            raise FileNotFoundError(f"Unable to locate the -r reads directory '{args.rDir}'")
        
    # Locate and validate that sequences exist in the sequences directory
    sequenceSuffix = ".cds" # always CDS for Salmon mapping
    args.sequenceFiles = validate_sequence_files(sequencesDir, sequenceSuffix)
    
    # Validate that salmon executable is locateable
    if args.salmonExe == None:
        salmon = distutils.spawn.find_executable("salmon")
        if salmon == None:
            raise FileNotFoundError("--salmon wasn't specified, and 'salmon' is missing from your PATH")
        args.salmonExe = salmon
    else:
        if not os.path.isfile(args.salmonExe):
            raise FileNotFoundError(f"Unable to locate the salmon executable '{args.salmonExe}'")
    args.salmonExe = os.path.abspath(args.salmonExe)
    
    # Validate numeric parameters
    if args.threads < 1:
        raise ValueError("--threads should be given a value >= 1")
    
    return locations

def validate_filter_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Derive the analysis directory name (without following symlink yet)
    locations.runName = args.analysisFolder
    
    # Validate that analysis directory exists
    args.runDirName = locations.resolve_runName(locations.analysisDir, locations.runName)
    
    # Validate that the cluster file exists
    resultsDir = os.path.join(analysisDir, args.runDirName)
    clusterFile = os.path.join(resultsDir, "BINge_clustering_result.tsv")
    if not os.path.isfile(clusterFile) and os.path.exists(clusterFile + ".ok"):
        raise FileNotFoundError("Unable to locate 'BINge_clustering_result.tsv' and " + 
                                f"'BINge_clustering_result.tsv.ok' within '{resultsDir}'")
    args.bingeFile = clusterFile
    
    # Validate cluster file format
    isBinge = validate_cluster_file(clusterFile)
    if not isBinge:
        raise ValueError(f"The file '{clusterFile}' does not appear to be a BINge cluster file; " + 
                         "has this been corrupted somehow?")
    
    # Validate BLAST file (if applicable)
    if args.useBLAST:
        blastDir = os.path.join(args.workingDirectory, "blast")
        blastFile = os.path.join(blastDir, f"MMseqs2_results.tsv")
        if not os.path.isfile(blastFile) and os.path.exists(blastFile + ".ok"):
            raise FileNotFoundError(f"Unable to locate 'MMseqs2_results.tsv' and " + 
                                    f"'MMseqs2_results.tsv.ok' within '{blastDir}'")
        args.blastFile = blastFile
    
    # Validate GFF3 file (if applicable)
    if args.useGFF3:
        sequencesDir = os.path.join(args.workingDirectory, "sequences")
        gff3Dir = os.path.join(sequencesDir, "gff3s")
        gff3Files = []
        for file in os.listdir(gff3Dir):
            if file.endswith(".gff3"):
                gff3Files.append(os.path.join(gff3Dir, file))
        if len(gff3Files) == 0:
            raise FileNotFoundError(f"Unable to locate any GFF3 files within '{gff3Dir}'")
        args.gff3Files = gff3Files
    
    # Validate salmon files (if applicable)
    if args.useSalmon:
        if args.salmonFiles != []:
            args.salmonFileFormat = validate_salmon_files(args.salmonFiles)
    
    # Validate numeric parameters
    if args.evalue < 0:
        raise ValueError("--evalue must be a positive float value >= 0")
    if args.readLength < 1:
        raise ValueError("--readLength must be a positive integer value >= 1")
    if args.minimumLength < 0:
        raise ValueError("--length must be a positive integer value >= zero")
    
    return locations

def validate_representatives_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Derive the analysis directory name (without following symlink yet)
    locations.runName = args.analysisFolder
    
    # Locate which analysis directory to use
    args.bingeFile = _locate_raw_or_filtered_results(args)
    
    # Validate cluster file format
    isBinge = validate_cluster_file(clusterFile)
    if not isBinge:
        raise ValueError(f"The file '{clusterFile}' does not appear to be a BINge cluster file; " + 
                         "has this been corrupted somehow?")
    
    # Validate numeric parameters
    if args.evalue < 0:
        raise ValueError("--evalue must be a positive float value >= 0")
    
    # Validate BLAST file (if applicable)
    if args.useBLAST:
        blastDir = os.path.join(args.workingDirectory, "blast")
        blastFile = os.path.join(blastDir, f"MMseqs2_results.tsv")
        if not os.path.isfile(blastFile) and os.path.exists(blastFile + ".ok"):
            raise FileNotFoundError(f"Unable to locate 'MMseqs2_results.tsv' and " + 
                                    f"'MMseqs2_results.tsv.ok' within '{blastDir}'")
        args.blastFile = blastFile
    
    # Validate GFF3 file (if applicable)
    if args.useGFF3:
        sequencesDir = os.path.join(args.workingDirectory, "sequences")
        gff3Dir = os.path.join(sequencesDir, "gff3s")
        gff3Files = []
        for file in os.listdir(gff3Dir):
            if file.endswith(".gff3"):
                gff3Files.append(os.path.join(gff3Dir, file))
        if len(gff3Files) == 0:
            raise FileNotFoundError(f"Unable to locate any GFF3 files within '{gff3Dir}'")
        args.gff3Files = gff3Files
    
    # Validate salmon files if relevant
    ## TBD: Implement properly
    #if args.salmonFiles != []:
    #    args.salmonFileFormat = validate_salmon_files(args.salmonFiles)
    
    return locations

def validate_dge_args(args):
    # Validate working directory
    locations = Locations(args.workingDirectory)
    
    # Derive the analysis directory name (without following symlink yet)
    locations.runName = args.analysisFolder
    
    # Locate which analysis directory to use
    args.bingeFile = _locate_raw_or_filtered_results(args)
    
    # Validate cluster file format
    isBinge = validate_cluster_file(clusterFile)
    if not isBinge:
        raise ValueError(f"The file '{clusterFile}' does not appear to be a BINge cluster file; " + 
                         "has this been corrupted somehow?")
    
    return locations

def _locate_raw_or_filtered_results(args, locations):
    rawRunDir = os.path.join(locations.analysisDir, locations.runName)
    filterRunDir = os.path.join(locations.filterDir, locations.runName)
    
    if os.path.isdir(filterRunDir):
        print(f"Detected '{args.runName}' within the 'filter' directory; will attempt to use this.")
        
        # Update the folder location (following symlink if necessary)
        filterRunDir = locations.resolve_runName(locations.filterDir)
        
        # Validate that the cluster file exists
        clusterFile = os.path.join(filterRunDir, locations.filteredClusterFile)
        if not os.path.isfile(clusterFile) and not os.path.exists(clusterFile + ".ok"):
            raise FileNotFoundError(f"Unable to locate '{locations.filteredClusterFile}' and " +
                                    f"'{locations.filteredClusterFile}.ok' within '{filterRunDir}'")
    else:
        print("No filtered results found; will try to use raw results instead.")
        print(f"Detected '{args.runName}' within the 'analysis' directory; will attempt to use this.")
        
        # Update the folder location (following symlink if necessary)
        rawRunDir = locations.resolve_runName(locations.analysisDir)
        
        # Validate that the cluster file exists
        clusterFile = os.path.join(rawRunDir, locations.clusterFile)
        if not os.path.isfile(clusterFile) and not os.path.exists(clusterFile + ".ok"):
            raise FileNotFoundError(f"Unable to locate '{locations.clusterFile}' and " +
                                    f"'{locations.clusterFile}.ok' within '{rawRunDir}'")
    return clusterFile

# File validations
def validate_salmon_files(salmonFiles):
    '''
    Validates Salmon files for 1) their existence and 2) their consistency of file format.
    Quits program if validation fails, so be warned!
    
    Parameters:
        salmonFiles -- a list containing strings pointing to Salmon files.
    Returns:
        fileFormat -- a string equal to "ec" if input files are equivalence classes,
                      or "quant" if they are quant.sf files.
    '''
    # Validate that salmon files exist
    for salmonFile in salmonFiles:
        if not os.path.isfile(salmonFile):
            raise FileNotFoundError(f"Unable to locate the salmon input file '{salmonFile}'")
    
    # Validate that input files are all of a consistent format
    isEC = False
    isQuant = False
    for salmonFile in salmonFiles:
        thisFileValid = False
        with open(salmonFile, "r") as fileIn:
            # Get the first 3 lines out of the file
            firstLine = fileIn.readline().rstrip("\r\n ")
            secondLine = fileIn.readline().rstrip("\r\n ")
            thirdLine = fileIn.readline().rstrip("\r\n ")
            
            # Check if it conforms to equivalence class expectations
            if firstLine.isdigit() and secondLine.isdigit() and not thirdLine.isdigit():
                isEC = True
                thisFileValid = True
            
            # Check if it conforms to quant file expectations
            elif firstLine.split("\t") == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
                isQuant = True
                thisFileValid = True
        
        if not thisFileValid:
            raise ValueError(f"The input file '{salmonFile}' does not appear to be a Salmon quant or " + 
                             "equivalence class file")
    
    if isEC and isQuant:
        raise ValueError("You appear to have given a mix of quant and equivalence class files; " +
                         "please only give one type.")
    
    return "ec" if isEC else "quant"

def validate_cluster_file(clusterFile):
    with open(clusterFile, "r") as fileIn:
        # Get the first two lines out of the file
        firstLine = fileIn.readline().rstrip("\r\n ")
        secondLine = fileIn.readline().rstrip("\r\n ")
        
        # Check if it conforms to BINge cluster file expectations
        if firstLine.startswith("#BINge clustering information file") \
            and secondLine.startswith("\t".join(["cluster_num", "sequence_id", "cluster_type"])):
                isBinge = True
        
        # Check if it conforms to CD-HIT file expectations
        elif firstLine.startswith(">Cluster 0") and secondLine.startswith("0"):
            isBinge = False
        
        # Raise an error otherwise
        else:
            errorMsg = (f"The input file '{clusterFile}' does not appear to be a BINge or " + 
                        "CD-HIT cluster file.\nYou should check your inputs and try again.")
            raise ValueError(errorMsg)

    return isBinge

def validate_fasta(fastaFile):
    '''
    Simple validator which just checks that the first character in a file is a '>'
    which likely means it's a FASTA file. If there are errors in the format that will
    not be detected here.
    
    Parameters:
        fastaFile -- a string indicating the location of a file to check.
    Returns:
        isFASTA -- a boolean where True means it is (probably) a FASTA, and False otherwise.
    '''
    with open(fastaFile, "r") as fileIn:
        firstLine = fileIn.readline()
        if not firstLine.startswith(">"):
            return False
        return True

def validate_sequence_files(sequencesDir, sequenceSuffix):
    sequenceFiles = [
        os.path.join(sequencesDir, f)
        for f in os.listdir(sequencesDir)
        if f.endswith(sequenceSuffix)
    ]
    if sequenceFiles == []:
        raise FileNotFoundError(f"Unable to locate any '{sequenceSuffix}' files within '{sequencesDir}'; " +
                                "have you run the initialisation step yet?")
    else:
        for sequenceFile in sequenceFiles:
            if not os.path.isfile(sequenceFile + ".ok"):
                raise FileNotFoundError(f"Unable to locate '{sequenceFile}.ok' to go along with '{sequenceFile}'; " +
                                        "re-run the initialisation step to ensure file is OK.")
    return sequenceFiles

def handle_symlink_change(existingLink, newLinkLocation):
    '''
    Detects if a symlink already exists and if it does, checks if it's pointing to the
    same location as the new link location. If it's not, an error is raised.
    
    Parameters:
        existingLink -- a string indicating the existing symlink file.
        newLink -- a string indicating where the symlink should point to.
    '''
    existingResolved = str(Path(existingLink).resolve())
    newResolved = str(Path(newLinkLocation).resolve())
    
    if existingResolved != newLinkLocation and existingResolved != newResolved:
        msg = f"File '{existingLink}' already points to '{existingResolved}'; and " + \
              f"you're trying to change it to '{newLinkLocation}'; initialise a new directory instead."
        raise FileExistsError(msg)

def touch_ok(fileName):
    '''
    Creates a file with the '.ok' suffix to indicate that a process has completed.
    
    Parameters:
        fileName -- a string indicating the file to create with the '.ok' suffix.
    '''
    if not fileName.endswith(".ok"):
        fileName = f"{fileName}.ok"
    open(fileName, "w").close()
