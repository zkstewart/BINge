import os
from pathlib import Path

class Locations:
    def __init__(self, workingDirectory):
        self.workingDirectory = workingDirectory
        
        # Set defaults
        self.runName = "most_recent"
    
    @property
    def workingDirectory(self):
        return self._workingDirectory
    
    @workingDirectory.setter
    def workingDirectory(self, value):
        value = os.path.abspath(value)
        if not os.path.isdir(value):
            raise FileNotFoundError(f"Unable to locate the working directory '{value}'")
        self._workingDirectory = value
    
    # Naive directory properties
    @property
    def analysisDir(self):
        return os.path.join(self.workingDirectory, "analysis")
    
    @property
    def tmpDir(self):
        return os.path.join(self.analysisDir, "tmp")
    
    @property
    def genomesDir(self):
        return os.path.join(self.workingDirectory, "genomes")
    
    @property
    def mappingDir(self):
        return os.path.join(self.workingDirectory, "mapping")
    
    @property
    def sequencesDir(self):
        return os.path.join(self.workingDirectory, "sequences")
    
    @property
    def txDir(self):
        return os.path.join(self.sequencesDir, "transcripts")
    
    @property
    def gff3Dir(self):
        return os.path.join(self.sequencesDir, "gff3s")
    
    @property
    def filterDir(self):
        return os.path.join(self.workingDirectory, "filter")
    
    @property
    def salmonDir(self):
        return os.path.join(self.workingDirectory, "salmon")
    
    @property
    def blastDir(self):
        return os.path.join(self.workingDirectory, "blast")
    
    @property
    def representativesDir(self):
        return os.path.join(self.workingDirectory, "representatives")
    
    @property
    def dgeDir(self):
        return os.path.join(self.workingDirectory, "dge")
    
    @property
    def annotateDir(self):
        return os.path.join(self.workingDirectory, "annotate")
    
    # Naive file properties
    @property
    def inputsJson(self):
        return "inputs.json"
    
    @property
    def clusterFile(self):
        return "BINge_clustering_result.tsv"
    
    @property
    def filteredClusterFile(self):
        return "BINge_clustering_result.filtered.tsv"
    
    @property
    def rScriptFile(self):
        return "BINge_DGE.R"
    
    @property
    def parametersFile(self):
        return "parameters.json"
    
    @property
    def pickleFile(self):
        return ".binge.pkl"
    
    @property
    def targetFile(self):
        return "targetFile.fasta"
    
    @property
    def blastFile(self):
        return "MMseqs2_results.tsv"
    
    @property
    def tx2geneFile(self):
        return "tx2gene.tsv"
    
    @property
    def salmonQCFile(self):
        return "salmon_qc_statistics.tsv"
    
    @property
    def salmonSampleFile(self):
        return "salmon_samples.txt"
    
    @property
    def annotationFile(self):
        return "BINge_annotation.tsv"
    
    @property
    def representativeMRNA(self):
        return "BINge_clustering_representatives.mrna"
    
    @property
    def representativeCDS(self):
        return "BINge_clustering_representatives.cds"
    
    @property
    def representativeAA(self):
        return "BINge_clustering_representatives.aa"
    
    # Smart properties
    @property
    def gff3Files(self, targetGenomes, annotatedGenomes):
        gff3Files = []
        
        # Check targetGenomes
        for targetGenome in targetGenomes:
            if targetGenome.gff3 != None:
                gff3Files.append(targetGenome.gff3)
        
        # Check annotatedGenomes
        for annotatedGenome in annotatedGenomes:
            gff3Files.append(annotatedGenome.gff3) # should always be non-null
        
        # There is no need to check that gff3Files != [], as it is possible to run BINge with no GFF3 inputs whatsoever
        return gff3Files
    
    @property
    def salmonFiles(self):
        # Locate salmon directories
        salmonDirs = [
            os.path.join(self.salmonDir, f)
            for f in os.listdir(self.salmonDir)
            if not f.endswith(".ok") # ignore any .ok files
            and not f.startswith("concatenated") # ignore the concatenated CDS FASTA and salmonDB
        ]
        
        # Validate that all directories have an associated .ok file
        for salmonDir in salmonDirs:
            if not os.path.exists(salmonDir + ".ok"):
                raise FileNotFoundError(f"Failed to find a .ok file associated with '{salmonDir}'")
        
        # Validate that salmon directories exist
        if len(salmonDirs) == 0:
            raise FileNotFoundError(f"Unable to locate any OK Salmon directories within '{self.salmonDir}'; " +
                                    "have you run the 'salmon' step yet?")
        
        # Locate and validate the salmon files within each directory
        salmonFiles = []
        for sDir in salmonDirs:
            # Validate file existence
            salmonFile = os.path.join(sDir, "quant.sf")
            if not os.path.isfile(salmonFile):
                raise FileNotFoundError(f"Unable to locate a 'quant.sf' file within '{sDir}'")
            
            # Validate file format
            isQuant = False
            with open(salmonFile, "r") as fileIn:
                firstLine = fileIn.readline().rstrip("\r\n ")
                if firstLine.split("\t") == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
                    isQuant = True
            if not isQuant:
                raise ValueError(f"The input file '{salmonFile}' does not appear to be a Salmon quant file")
            
            # Store the file
            salmonFiles.append(salmonFile)
        
        return salmonFiles
    
    @property
    def runName(self):
        return self._runName
    
    @runName.setter
    def runName(self, value):
        assert isinstance(value, str)
        
        if value != "most_recent":
            if value.startswith("run_"):
                value = value
            else:
                value = f"run_{value}"
        else:
            value = value
        
        self._runName = value
    
    def resolve_runName(self, value):
        value = os.path.abspath(value)
        if not os.path.isdir(value):
            raise FileNotFoundError(f"Unable to locate '{value}'")
        
        runDir = os.path.join(value, self.runName)
        if not os.path.isdir(runDir):
            raise FileNotFoundError(f"Unable to locate '{self.runName}' within '{value}'")
        
        if self.runName != "most_recent":
            return runDir
        else:
            resolvedDir = str(Path(runDir).resolve()) # go from 'most_recent' to 'run_XX'
            if not os.path.isdir(resolvedDir):
                raise FileNotFoundError(f"Unable to locate '{os.path.basename(resolvedDir)}' within '{value}'")
            
            print(f"Run folder identified as: 'most_recent' -> '{resolvedDir}'")
            return resolvedDir
    
    def get_sequenceFiles(self, targetGenomes, annotatedGenomes, transcriptomes, sequenceSuffix):
        ACCEPTED_SUFFIXES = ["mrna", "cds", "aa"]
        if not sequenceSuffix in ACCEPTED_SUFFIXES:
            raise ValueError(f"get_sequenceFiles should receive a value in '{ACCEPTED_SUFFIXES}', not '{sequenceSuffix}'")
        
        sequenceFiles = []
        for inputList in [targetGenomes, annotatedGenomes, transcriptomes]:
            for inputObj in inputList:
                if sequenceSuffix == "mrna":
                    if inputObj.mrna != None:
                        sequenceFiles.append(inputObj.mrna)
                elif sequenceSuffix == "cds":
                    if inputObj.cds != None:
                        sequenceFiles.append(inputObj.cds)
                else:
                    if inputObj.aa != None:
                        sequenceFiles.append(inputObj.aa)
        
        if sequenceFiles == []:
            raise FileNotFoundError(f"Unable to locate any '{sequenceSuffix}' files; " +
                                    f"have you run the initialisation step yet?")
        else:
            for sequenceFile in sequenceFiles:
                if not os.path.isfile(sequenceFile + ".ok"):
                    raise FileNotFoundError(f"Unable to locate '{sequenceFile}.ok' to go along with '{sequenceFile}'; " +
                                            "re-run the initialisation step to ensure file is OK.")
        return sequenceFiles
