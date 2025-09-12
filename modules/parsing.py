import os, re, pickle, codecs, gzip
from contextlib import contextmanager

from .salmon import EquivalenceClassCollection, QuantCollection, DGEQuantCollection

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class BLAST_Results:
    '''
    The BLAST_Results Class provides parsing capability for outfmt6 files.
    
    Attributes:
        file -- a string indicating the location of the outfmt6 file that was
                parsed to generate this object instance
        evalue -- an int or float value indicating the cutoff to enforce for
                  retrieving hits for a query
        num_hits -- an int value which can be -1 (to retrieve all hits) or a value
                    > 0 to get that many hits per query
    '''
    def __init__(self, outfmt6File):
        # Validate input types
        assert isinstance(outfmt6File, str)
        
        # Validate that input file location exists
        if not os.path.isfile(outfmt6File):
            raise Exception(f"'{outfmt6File} does not point to an existing file location")
        self.file = outfmt6File
        
        # Set default property values
        self.evalue = 1e-5
        self.num_hits = -1
        self.results = None
        
        # Also set helper attribute
        self.isBLAST_Results = True
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int), \
            "evalue can only be set to a float or integer value"
        
        if value < 0:
            raise ValueError("evalue must be a value >= 0")
        
        self._evalue = value
    
    @property
    def num_hits(self):
        return self._num_hits
    
    @num_hits.setter
    def num_hits(self, value):
        assert isinstance(value, int), \
            "num_hits must be an integer value"
        
        if value == 0:
            raise ValueError("num_hits cannot be 0 (which would retrieve no hits)")
        if value < -1:
            raise ValueError("num_hits cannot be less than -1 (which retrieves all hits)")
        
        self._num_hits = value
    
    def parse(self):
        '''
        Parameters:
            self.file -- a string indicating the location of a BLAST results file in
                         outfmt6 format.
        Sets:
            self.results -- a dict with structure like:
                                query_id1: [
                                    [target_id, identity_pct, query_start, query_end, target_start, target_end, evalue],
                                    ...
                                ],
                                query_id2: [ ... ],
                                ...
        '''
        blastDict = {}
        
        with open(self.file, "r", encoding=get_codec(self.file)) as fileIn:
            for line in fileIn:
                # Extract details
                sl = line.rstrip("\r\n").split("\t")
                qid = sl[0]
                tid = sl[1]
                identityPct = float(sl[2])
                qstart = int(sl[6])
                qend = int(sl[7])
                tstart = int(sl[8])
                tend = int(sl[9])
                evalue = float(sl[10])
                bitscore = float(sl[11])
                
                # Skip if evalue isn't significant
                if evalue > self.evalue: # filter here since self.evalue might differ between BLAST run and parsing now
                    continue
                
                # Store result
                if qid not in blastDict:
                    blastDict[qid] = [[tid, identityPct, qstart, qend, tstart, tend, evalue, bitscore]]
                else:
                    blastDict[qid].append([tid, identityPct, qstart, qend, tstart, tend, evalue, bitscore])
        
        # Sort individual entries in blastDict
        for value in blastDict.values():
            value.sort(key = lambda x: (x[6], -x[7])) # sort by evalue (lower) and bitscore (higher)
        
        # Enforce num_hits threshold
        if self.num_hits != -1:
            for key in blastDict.keys():
                blastDict[key] = blastDict[key][0:self.num_hits]
        
        self.results = blastDict
    
    def __iter__(self):
        return iter(self.results)
    
    def __len__(self):
        return len(self.results)
    
    def __getitem__(self, key):
        if key in self.results:
            return self.results[key]
    
    def __contains__(self, item):
        return True if self[item] is not None else False
    
    def __str__(self):
        return "BLAST_Results; parsed '{0}' file, contains results for {1} queries".format(
            self.file, len(self)
        )
    
    def __repr__(self):
        return "BLAST_Results(outfmt6='{0}', queries={1}, evalue={2}, num_hits={3})".format(
            self.file, len(self), self.evalue, self.num_hits
        )

class BINge_Results:
    seqFileRegex = re.compile(r"#sequence_file_(\d+)=(.+)")
    
    def __init__(self):
        self._binned = {}
        self._unbinned = {}
        self._highestBinnedNumber = -1
        self._highestUnbinnedNumber = -1
    
    @property
    def binned(self):
        return self._binned
    
    @binned.setter
    def binned(self, value):
        if not isinstance(value, dict):
            raise TypeError(f".binned must be a dict, not '{type(value)}'")
        
        if set(value.keys()).intersection(set(self.unbinned.keys())) != set():
            raise ValueError("Cluster ID collision between new .binned and existing .unbinned values!")
        
        highestNumber = -1
        for k, v in value.items():
            if not isinstance(k, int):
                raise TypeError(f".binned keys must be integer, not '{type(k)}'")
            highestNumber = max(highestNumber, k)
        
        self._binned = value
        self._highestBinnedNumber = highestNumber
    
    @property
    def unbinned(self):
        return self._unbinned
    
    @unbinned.setter
    def unbinned(self, value):
        if not isinstance(value, dict):
            raise TypeError(f".unbinned must be a dict, not '{type(value)}'")
        
        if set(value.keys()).intersection(set(self.binned.keys())) != set():
            raise ValueError("Cluster ID collision between new .unbinned and existing .binned values!")
        
        highestNumber = -1
        for k, v in value.items():
            if not isinstance(k, int):
                raise TypeError(f".unbinned keys must be integer, not '{type(k)}'")
            highestNumber = max(highestNumber, k)
        
        self._unbinned = value
        self._highestUnbinnedNumber = highestNumber
    
    def update_unbinned_ids(self, unbinnedClusterDict):
        '''
        Receives a cluster dictionary as obtained by parsing an external clusterer's results
        e.g., MM_Clust.parse_tsv, where cluster keys will start at index 0. To prevent key
        collision, you should run it through this function before attempting to set
        it as the .unbinned value.
        
        Parameters:
            unbinnedClusterDict -- a dictionary where keys are integers and values are
                                   lists of sequence IDs (strings)
        '''
        return { k+self._highestBinnedNumber+1:v for k,v in unbinnedClusterDict.items() }
    
    def parse(self, bingeFile):
        '''
        Reads in the output file of BINge as a dictionary assocating clusters to their
        sequence members.
        
        Parameters:
            bingeFile -- a string pointing to the location of a BINge cluster output file.
        Sets:
            self.binned / self.unbinned -- a dictionary with structure like:
                        {
                            0: [seqid1, seqid2, ...],
                            1: [ ... ],
                            ...
                        }
        '''
        self.binned = {}
        self.unbinned = {}
        firstCommented = True
        firstUncommented = True
        with open(bingeFile, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t ")
                
                # Handle header lines
                if line.startswith("#"):
                    # Handle file format identifier line
                    if firstCommented:
                        assert line.startswith("#BINge clustering information file"), \
                            ("BINge file is expected to start with a specific comment line! " +
                            "Your file is hence not recognised as a valid BINge cluster file.")
                        firstCommented = False
                    # Handle sequence file identifier lines
                    else:
                        continue # silently skip any other comments a user might have added in
                # Handle blank lines
                elif l == "": # try to allow whitespace lines if they've somehow been added in
                    continue
                # Handle body lines
                else:
                    sl = line.rstrip("\r\n ").split("\t")
                    
                    # Handle column labels
                    if firstUncommented:
                        assert sl == ["cluster_num", "sequence_id", "cluster_type"], \
                            ("BINge file is expected to have a specific header line on the first uncommented line! " +
                            "Your file is hence not recognised as a valid BINge cluster file.")
                        firstUncommented = False
                    # Handle all subsequent information-containing lines
                    else:
                        clustNum, seqID, clusterType = int(sl[0]), sl[1], sl[2]
                        
                        if clusterType == "binned":
                            self.binned.setdefault(clustNum, [])
                            self.binned[clustNum].append(seqID)
                            if clustNum > self._highestBinnedNumber:
                                self._highestBinnedNumber = clustNum
                        elif clusterType == "unbinned":
                            self.unbinned.setdefault(clustNum, [])
                            self.unbinned[clustNum].append(seqID)
                            if clustNum > self._highestUnbinnedNumber:
                                self._highestUnbinnedNumber = clustNum
                        else:
                            raise ValueError(f"Did not recognise '{clusterType}' as a valid third column value in a BINge results file")
    
    def filter(self, clustersToRemove):
        '''
        Drops any clusters from .binned or .unbinned with a key value match to
        clusters listed in clustersToRemove. Providing a key not found
        in this object will raise an error.
        
        Parameters:
            clustersToRemove -- a list or set of cluster IDs to drop from this object's
                                .binned or .unbinned dictionaries. 
        '''
        # Locate and validate all IDs for removal
        binnedToRemove = []
        unbinnedToRemove = []
        for clusterID in clustersToRemove:
            if clusterID in self.binned:
                binnedToRemove.append(clusterID)
            elif clusterID in self.unbinned:
                unbinnedToRemove.append(clusterID)
            else:
                raise ValueError(f"No cluster with number {clusterID} found in this BINge results object")
        
        # Convert to set for lookup
        binnedToRemove = set(binnedToRemove)
        unbinnedToRemove = set(unbinnedToRemove)
        
        # Drop clusters
        self.binned = { k:v for k,v in self.binned.items() if not k in binnedToRemove }
        self.unbinned = { k:v for k,v in self.unbinned.items() if not k in unbinnedToRemove }
    
    def write(self, outputFileName, clusterTypes="all"):
        '''
        Generates a BINge results file based on the values in self.binned and self.unbinned
        
        Parameters:
            outputFileName -- a string indicating the location to write the output file
            clusterTypes -- a string in the list of ["all", "binned", "unbinned"] signalling which
                            result type(s) to write to output
        '''
        ACCEPTED_TYPES = ["all", "binned", "unbinned"]
        if not clusterTypes in ACCEPTED_TYPES:
            raise ValueError(f"write() does not recognise '{clusterTypes}' as a valid value for clusterTypes")
        
        with open(outputFileName, "w") as fileOut:
            # Write header lines
            fileOut.write("#BINge clustering information file\n")
            fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
            
            # Write content lines
            if (clusterTypes == "all" or clusterTypes == "binned") and hasattr(self, "binned"):
                for clusterNum, seqValues in self.binned.items():
                    for seqID in seqValues:
                        fileOut.write(f"{clusterNum}\t{seqID}\tbinned\n")
            if (clusterTypes == "all" or clusterTypes == "unbinned") and hasattr(self, "unbinned"):
                for clusterNum, seqValues in self.unbinned.items():
                    for seqID in seqValues:
                        fileOut.write(f"{clusterNum}\t{seqID}\tunbinned\n")
    
    def __iter__(self):
        for key, value in self.binned.items():
            yield key, value
        for key, value in self.unbinned.items():
            yield key, value
    
    def __len__(self):
        return len(self.binned) + len(self.unbinned)
    
    def __repr__(self):
        return "<BINge_Results object;num_binned_clusters={0};num_unbinned_clusters={1};total_clustered_sequences={2}>".format(
            len(self.binned),
            len(self.unbinned),
            sum([ len(values) for values in self.binned.values() ]) + sum([ len(values) for values in self.unbinned.values() ])
        )

def parse_equivalence_classes(equivalenceClassFiles, sampleNames):
    '''
    Parses in one or more equivalence class files from Salmon, producing
    an EquivalenceClassCollection object enabling read count summarisation.
    
    Parameters:
        equivalenceClassFiles -- a list containing strings pointing to the location
                                 of Salmon equivalence class files
                                 (eq_classes.txt files).
        sampleNames -- an equal length list indicating the sample names for each
                       equivalence class file.
    '''
    assert len(equivalenceClassFiles) == len(sampleNames), \
        ("parse_equivalence_classes cannot parse equivalence classes since the number " +
         f"of files ({len(equivalenceClassFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    ecCollection = EquivalenceClassCollection()
    for i in range(len(equivalenceClassFiles)):
        eqFile = equivalenceClassFiles[i]
        sample = sampleNames[i]
        
        ecCollection.parse_eq_file(eqFile, sample)
    return ecCollection

def parse_quants(quantFiles, sampleNames):
    '''
    Parses in one or more quant files from Salmon, producing a QuantCollection
    object enabling read count summarisation.
    
    Parameters:
        quantFiles -- a list containing strings pointing to the location
                      of Salmon quant files (quant.sf files).
        sampleNames -- an equal length list indicating the sample names for each
                       salmon quant file.
    '''
    assert len(quantFiles) == len(sampleNames), \
        ("parse_quants cannot parse quant files since the number " +
         f"of files ({len(quantFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    quantCollection = QuantCollection()
    for i in range(len(quantFiles)):
        quantFile = quantFiles[i]
        sample = sampleNames[i]
        
        quantCollection.parse_quant_file(quantFile, sample)
    return quantCollection

def parse_dge_quants(quantFiles, sampleNames):
    '''
    Parses in one or more quant files from Salmon, producing a DGEQuantCollection
    object enabling the generation of output files needed for DESeq2 DGE analysis.
    
    Parameters:
        quantFiles -- a list containing strings pointing to the location
                      of Salmon quant files (quant.sf files).
        sampleNames -- an equal length list indicating the sample names for each
                       salmon quant file.
    '''
    assert len(quantFiles) == len(sampleNames), \
        ("parse_dge_quants cannot parse quant files since the number " +
         f"of files ({len(quantFiles)}) does not match the number of " +
         f"sample names ({len(sampleNames)})")
    
    dgeQuantCollection = DGEQuantCollection()
    for i in range(len(quantFiles)):
        quantFile = quantFiles[i]
        sample = sampleNames[i]
        
        dgeQuantCollection.parse_quant_file(quantFile, sample)
    return dgeQuantCollection

def parse_gff3_ids(gff3File):
    '''
    Simple function to iterate through a GFF3 file and hold onto any parent and
    subfeature IDs. This process should encompass genes and mRNAs as well as any
    other subfeature e.g., lnc_RNA, whilst skipping over exons and CDS'.
    
    Parameters:
        gff3File -- a string indicating the location of a GFF3 file.
    Returns:
        annotIDs -- a set of strings containing the IDs of all parent and subfeatures
                    found in the GFF3 file.
    '''
    # Setup for slim parsing of attributes
    def _format_attributes(attributes):
        "Code borrowed from ZS_GFF3IO"
        SLIM_ATTRIBUTES = ["id", "parent"]
        
        splitAttributes = []
        for a in attributes.split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        attributesDict = {
            splitAttributes[i]: splitAttributes[i+1]
            for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)
            if splitAttributes[i].lower() in SLIM_ATTRIBUTES
        }
        return attributesDict
    
    parentIDs = set()
    subIDs = set()
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Parse content lines
            if not line.startswith("#") and len(sl) == 9:
                featureType, attributes = sl[2], sl[8]
                attributesDict = _format_attributes(attributes)
                if featureType == "gene":
                    parentIDs.add(attributesDict["ID"])
                elif "Parent" in attributesDict and attributesDict["Parent"] in parentIDs:
                    subIDs.add(attributesDict["ID"])
    annotIDs = parentIDs.union(subIDs)
    return annotIDs

def load_sequence_length_index(indexFile):
    '''
    Load in the pickled result of generate_sequence_length_index().
    
    Parameters:
        indexFile -- a string indicating the location of index generated from a FASTA file.
    '''
    if os.path.exists(indexFile):
        with open(indexFile, "rb") as fileIn:
            seqLenDict = pickle.load(fileIn)
        return seqLenDict
    else:
        raise FileNotFoundError(("load_sequence_length_index() failed because " + 
                                 f"'{indexFile}' doesn't exist!"))

def locate_read_files(readsDir, readsSuffix, isSingleEnd):
    '''
    Parameters:
        readsDir -- a list of strings indicating the directory where
                    reads files are located.
        readsSuffix -- a string indicating the suffix of the reads files.
        isSingleEnd -- a boolean indicating whether the reads are single-end (True)
                       or paired-end (False).
    Returns:
        forwardReads -- a list of strings indicating the location of forward reads.
        reverseReads -- a list of strings indicating the location of reverse reads
                        OR None if the reads are single-end.
        sampleNames -- a list of strings indicating the common prefix of the reads files
                       (if paired-end) or the file name sans suffix (if single-end).
    '''
    # Locate files from the directory
    forwardReads = []
    reverseReads = []
    sampleNames = []
    for dir in readsDir:
        for file in os.listdir(dir):
            if file.endswith(readsSuffix):
                if isSingleEnd:
                    forwardReads.append(os.path.join(dir, file))
                    sampleNames.append(file.rsplit(readsSuffix, maxsplit=1)[0])
                else:
                    if file.endswith(f"1{readsSuffix}"):
                        forwardReads.append(os.path.join(dir, file))
                    elif file.endswith(f"2{readsSuffix}"):
                        reverseReads.append(os.path.join(dir, file))
                    else:
                        raise ValueError(f"{file} ends with the expected suffix '{readsSuffix}' but is " +
                                        "not preceeded by a 1 or 2!")
    forwardReads.sort()
    reverseReads.sort()
    
    # Validate that paired files match
    if not isSingleEnd:
        if not len(forwardReads) == len(reverseReads):
            raise FileNotFoundError(f"Number of reads don't match for forward ({len(forwardReads)}) " + 
                                    f"and reverse ({len(reverseReads)}) files")
        
        for i in range(len(forwardReads)):
            forward = os.path.basename(forwardReads[i]) # strip directory which may falsely give a common prefix
            reverse = os.path.basename(reverseReads[i])
            
            prefix = os.path.commonprefix([forward, reverse]).rstrip("_ ")
            if prefix == "":
                raise ValueError("Forward and reverse read pairs don't have a common prefix after " +
                                 "sorting; make sure that files are named in a way that they can be " +
                                 "paired together!")
            if not forwardReads[i].startswith(f"{prefix}1") and reverseReads[i].startswith(f"{prefix}2"):
                raise ValueError(f"Forward read '{forwardReads[i]}' and reverse read '{forwardReads[i]}' " +
                                 "files don't start with a valid common prefix i.e., " + 
                                 f"'{prefix} should preceed a 1 or 2")
            sampleNames.append(prefix)
    
    # Return files
    return forwardReads, reverseReads if reverseReads != [] else None, sampleNames

def parse_binge_representatives(representativesFasta):
    '''
    Reads in the output file of BINge 'representatives' to associate cluster IDs to
    their representative sequence ID.
    
    Parameters:
        representativesFasta -- a string pointing to the location of a BINge representatives
                                FASTA file.
    Returns:
        clustToRep -- a dictionary with cluster IDs as keys and representative sequence IDs
                      as values.
        repToClust -- a dictionary with representative sequence IDs as keys and cluster IDs
                      as values.
    '''
    clustToRep = {}
    repToClust = {}
    with open(representativesFasta, "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                clusterID, repID = line[1:].rstrip("\r\n").split(" ", maxsplit=1)
                repID = repID.split("=", maxsplit=1)[1]
                
                clustToRep[clusterID] = repID
                repToClust[repID] = clusterID
    return clustToRep, repToClust
