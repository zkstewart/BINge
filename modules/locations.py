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
    def gff3Files(self):
        gff3Files = []
        
        # Check self.gff3Dir
        for file in os.listdir(self.gff3Dir):
            if file.endswith(".gff3"):
                gff3Files.append(os.path.join(self.gff3Dir, file))
        
        # Check self.genomesDir
        for file in os.listdir(self.genomesDir):
            if file.endswith(".gff3"):
                gff3Files.append(os.path.join(self.genomesDir, file))
        
        if len(gff3Files) == 0:
            raise FileNotFoundError(f"Unable to locate any GFF3 files within '{self.gff3Dir}' " + 
                                    f"or '{self.genomesDir}'")
        return gff3Files
    
    @property
    def salmonFiles(self):
        # Locate salmon directories
        salmonOK = [
            os.path.join(self.salmonDir, f)
            for f in os.listdir(self.salmonDir)
            if f.endswith(".ok")
            and not ".salmonDB" in f # ignore the salmonDB directory
            and not f.startswith("concatenated") # ignore the concatenated CDS file
        ]
        salmonDirs = [
            f.rsplit(".ok", maxsplit=1)[0]
            for f in salmonOK
        ]
        
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
    
    def get_sequenceFiles(self, sequenceSuffix):
        sequenceFiles = [
            os.path.join(self.sequencesDir, f)
            for f in os.listdir(self.sequencesDir)
            if f.endswith(sequenceSuffix)
        ] + [
            os.path.join(self.genomesDir, f)
            for f in os.listdir(self.genomesDir)
            if f.endswith(sequenceSuffix)
        ]
        if sequenceFiles == []:
            raise FileNotFoundError(f"Unable to locate any '{sequenceSuffix}' files within '{self.sequencesDir}' " +
                                    f"or '{self.sequencesDir}'; have you run the initialisation step yet?")
        else:
            for sequenceFile in sequenceFiles:
                if not os.path.isfile(sequenceFile + ".ok"):
                    raise FileNotFoundError(f"Unable to locate '{sequenceFile}.ok' to go along with '{sequenceFile}'; " +
                                            "re-run the initialisation step to ensure file is OK.")
        return sequenceFiles
