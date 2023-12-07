import os, distutils.spawn, platform, subprocess, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages.ZS_Utility import convert_wsl_to_windows_path

def find_executable(exeName):
    '''
    Searches for the given executable name in the system PATH. If found, it will
    return the full path to it. Otherwise, it will return None.
    
    Parameter:
        exeName -- a string of an executable file name.
    Returns:
        which -- a string similar to that returned by Linux's which, or None if not found.
    '''
    if platform.system() != "Windows":
        return distutils.spawn.find_executable(exeName)
    else:
        cmd = ["wsl", "~", "-e", "which", exeName]        
        run_find_exe = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        findout, finderr = run_find_exe.communicate()
        
        if findout.decode("utf-8").rstrip("\r\n ") != "":
            wslWhich = findout.decode("utf-8").rstrip("\r\n ")
            which = convert_wsl_to_windows_path(wslWhich)
            return which
        else:
            return None

def validate_args(args):
    # Validate input file locations
    for inputArgument in args.inputFiles:
        assert inputArgument.count(",") <= 1, \
            print(f"-i value '{inputArgument}' has more than one comma in it; format is unhandled")
        
        for inputFile in inputArgument.split(","):
            if not os.path.isfile(inputFile):
                print(f'I am unable to locate the -i input file ({inputFile})')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
    
    for inputArgument in args.genomeFiles:
        assert inputArgument.count(",") <= 1, \
            print(f"-g value '{inputArgument}' has more than one comma in it; format is unhandled")
        
        for inputFile in inputArgument.split(","):
            if not os.path.isfile(inputFile):
                print(f'I am unable to locate the -g input file ({inputFile})')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
    
    # Validate optional BINge parameters
    if args.threads < 1:
        print("--threads should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if args.convergenceIters < 1:
        print("--convergence_iters should be given a value >= 1")
        print("Fix this and try again.")
        quit()
    if not 0.0 < args.identity <= 1.0:
        print("--identity should be given a value greater than zero, and equal to " + 
              "or less than 1")
        print("Fix this and try again.")
        quit()
    
    # Specifically handle gmapIdentity
    if args.gmapIdentity == []:
        args.gmapIdentity = [0.95 for _ in range(len(args.inputFiles))]
    if len(args.gmapIdentity) != len(args.inputFiles):
        print("--gmapIdentity parameter must have the same number of values as given to -i")
        print("Fix this and try again.")
        quit()
    for gmapID in args.gmapIdentity:
        if not 0.0 < gmapID <= 1.0:
            print("--gmapIdentity should be given values greater than zero, and equal to " + 
                "or less than 1")
            print("Fix this and try again.")
            quit()
    
    # Validate GMAP location
    if args.gmapDir == None:
        gmap = find_executable("gmap")
        gmap_build = find_executable("gmap_build")
        
        # Raise error if we failed to find the executables
        if gmap == None or gmap_build == None:
            print("--gmapDir wasn't specified, so I tried looking for it in your PATH.")
            print("I failed to find 'gmap' and/or 'gmap_build' by doing so.")
            print("Make sure they're in your PATH or explicitly indicate them with --gmapDir")
            quit()
        else:
            gmapDir = os.path.dirname(gmap)
            gmap_buildDir = os.path.dirname(gmap_build)
            
            # Raise errors if they're not in the same spot
            if gmapDir != gmap_buildDir:
                print("--gmapDir wasn't specified, so I tried looking for it in your PATH.")
                print("I found 'gmap' and 'gmap_build', but they're in different locations.")
                print("It would make my life a lot easier if you had them in the same spot.")
                print("Please fix that and then try again.")
                quit()
            
            # Set the value so we can use it later
            args.gmapDir = gmapDir
    else:
        if not os.path.isdir(args.gmapDir):
            print(f'I am unable to locate the GMAP directory ({args.gmapDir})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
    # Validate optional MMseqs2 parameters
    if args.unbinnedClusterer in ["mmseqs-cascade", "mmseqs-linclust"]:
        if args.mmseqsDir == None:
            print(f"--mmseqs must be specified if you're using '{args.unbinnedClusterer}'")
            quit()
        if not os.path.isdir(args.mmseqsDir):
            print(f'I am unable to locate the MMseqs2 directory ({args.mmseqsDir})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        "--tmpDir is validated by the MM_DB Class"
        if args.evalue < 0:
            print("--evalue must be greater than or equal to 0")
            quit()
        if not 0 <= args.coverage <= 1.0:
            print("--coverage must be a float in the range 0.0 -> 1.0")
            quit()
        "--mode is controlled by argparse choices"
        
        # Validate "MMS-CASCADE" parameters
        if args.sensitivity in ["5.7", "7.5"]:
            args.sensitivity = float(args.sensitivity)
        else:
            args.sensitivity = int(args.sensitivity)
        
        if args.steps < 1:
            print("--steps must be greater than or equal to 1")
            quit()
    
    # Validate optional CD-HIT parameters
    if not 0.0 <= args.cdhitShortCov <= 1.0:
        print("--cdhit_shortcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not 0.0 <= args.cdhitLongCov <= 1.0:
        print("--cdhit_longcov should be given a value in the range of 0 -> 1 (inclusive)")
        print("Fix this and try again.")
        quit()
    if not args.cdhitMem >= 100:
        print("--cdhit_mem should be given a value at least greater than 100 (megabytes)")
        print("Fix this and try again.")
        quit()
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory):
        print(f'Output directory specified already exists ({args.outputDirectory})')
        print('This means BINge will try to resume a previously started run.')
        print('No files will be overwritten, so if something went wrong previously you ' + 
              'should fix that up before running BINge again.')
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
            print(f'I am unable to locate the salmon input file ({salmonFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    
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
            print(f"The input file '{salmonFile}' does not appear to be a Salmon quant or " + 
                  "equivalence class file")
            print("You should check your inputs and try again.")
            quit()
    
    if isEC and isQuant:
        print("You appear to have given a mix of quant and equivalence class files.")
        print("That's too hard for me to figure out, so please only give one type and try again.")
        quit()
    
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
