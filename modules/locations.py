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
            raise FileNotFoundError(f"Unable to locate the working directory '{args.workingDirectory}'")
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
    
    # Naive file properties
    @property
    def clusterFile(self):
        return "BINge_clustering_result.tsv"
    
    @property
    def filteredClusterFile(self):
        return "BINge_clustering_result.filtered.tsv"
    
    @property
    def parametersFile(self):
        return "parameters.json"
    
    @property
    def pickleFile(self):
        return ".binge.pkl"
    
    @property
    def blastFile(self):
        return "MMseqs2_results.tsv"
    
    # Smart properties
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
        if self.runName != "most_recent":
            return self.runName
        else:
            value = os.path.abspath(value)
            if not os.path.isdir(value):
                raise FileNotFoundError(f"Unable to locate '{value}'")
            
            runLink = os.path.join(value, self.runName)
            if not os.path.isdir(runLink):
                raise FileNotFoundError(f"Unable to locate '{self.runName}' within '{value}'")
            
            resolvedDir = str(Path(runLink).resolve()) # go from 'most_recent' to 'run_XX'
            if not os.path.isdir(resolvedDir):
                raise FileNotFoundError(f"Unable to locate '{os.path.basename(resolvedDir)}' within '{value}'")
            
            print(f"Run folder identified as: 'most_recent' -> '{resolvedDir}'")
            return resolvedDir
