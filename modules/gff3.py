#! python3

import os, re, sys
from collections import Counter

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file
from gff3fasta import Sequence

class GFF3Feature:
    IMMUTABLE = ["ID", "ftype"] # these attributes should never change once set
    HEADER_FORMAT = ["contig", "source", "ftype", "start", "end", "score", "strand", "frame", "attributes"]
    CONTROLLED_ATTRIBUTES = ["ID", "Parent"] # these attributes are controlled by the GFF3Feature class
    
    def __init__(self, ID, ftype, start=None, end=None, strand=None, contig=None,
                 source=None, score=None, frame=None, attributes=None, children=None, parents=None,
                 isInferred=False):
        self.ID = ID
        self.ftype = ftype
        self.start = start
        self.end = end
        self.strand = strand
        self.contig = contig
        self.source = source
        self.score = score
        self.frame = frame
        self.attributes = attributes
        
        self._children = []
        self.children = children
        self.parents = parents
        self.isInferred = isInferred # flag for checking if this feature was parsed or inferred
        self.isGFF3Feature = True # flag for easier type checking
    
    @staticmethod
    def make_ftype_case_appropriate(ftype):
        if ftype.lower() == "gene":
            return "gene"
        elif ftype.lower() == "mrna":
            return "mRNA"
        elif ftype.lower() == "exon":
            return "exon"
        elif ftype.lower() == "cds":
            return "CDS"
        elif ftype.lower() == "lnc_rna":
            return "lnc_RNA"
        elif ftype.lower() == "product":
            return "Product"
        else:
            return ftype
    
    @property
    def start(self):
        if self._start != None:
            return self._start
        else:
            if len(self.children) == 0:
                return None
            return min([ child.start for child in self.children ])
    
    @start.setter
    def start(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"Start value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"Start value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._start = value
        else:
            self._start = None
    
    @property
    def end(self):
        if self._end != None:
            return self._end
        else:
            if len(self.children) == 0:
                return None
            return max([ child.end for child in self.children ])
    
    @end.setter
    def end(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"End value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"End value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._end = value
        else:
            self._end = None
    
    @property
    def strand(self):
        if self._strand != None and self._strand != ".":
            return self._strand
        else:
            childStrands = [ child.strand for child in self.children if child.strand != None and child.strand != "." ]
            if len(childStrands) == 0:
                return "+" # default strand if no children have a strand
            else:
                mostCommonStrand = Counter(childStrands).most_common(1)[0][0] # basic majority vote
                return mostCommonStrand
    
    @strand.setter
    def strand(self, value):
        ACCEPTED_STRANDS = ["+", "-", "."] # "." might represent unknown strand
        if value != None:
            if not value in ACCEPTED_STRANDS:
                raise ValueError(f"Strand value '{value}' is not recognised; should be one of {ACCEPTED_STRANDS}")
            self._strand = value
        else:
            self._strand = None
    
    @property
    def source(self):
        if self._source != None:
            return self._source
        else:
            return "."
    
    @source.setter
    def source(self, value):
        if value == None:
            pass
        elif not isinstance(value, str):
            raise ValueError(f"GFF3 source must be a string, not {type(value)}")
        self._source = value
    
    @property
    def score(self):
        if self._score != None:
            return self._score
        else:
            return "."
    
    @score.setter
    def score(self, value):
        self._score = value
    
    @property
    def frame(self):
        if self._frame != None:
            return self._frame
        else:
            return "."
    
    @frame.setter
    def frame(self, value):
        ACCEPTED_FRAMES = ["0", "1", "2", "."] # "." might represent unknown or irrelevant frame
        if value != None:
            if not value in ACCEPTED_FRAMES:
                raise ValueError(f"Frame value '{value}' is not recognised; should be one of {ACCEPTED_FRAMES}")
            self._frame = value
        else:
            self._frame = None
    
    @property
    def attributes(self):
        attrDict = {}
        
        # Set controlled attributes
        attrDict["ID"] = self.ID
        if len(self.parents) != 0:
            attrDict["Parent"] = ",".join(self.parents)
        
        # Add other attributes
        if self._attributes != None:
            for key, value in self._attributes.items():
                if not key in GFF3Feature.CONTROLLED_ATTRIBUTES:
                    attrDict[key] = value
        
        return attrDict
    
    @attributes.setter
    def attributes(self, value):
        if value == None:
            self._attributes = None
        else:
            if not isinstance(value, dict):
                raise TypeError(f"Attributes value must a dict, not '{type(value)}'")
            cleanedValue = { k:v for k,v in value.items() if not k in GFF3Feature.CONTROLLED_ATTRIBUTES }
            
            self._attributes = cleanedValue
    
    @property
    def children(self):
        return self._children
    
    @children.setter
    def children(self, value):
        if isinstance(value, list):
            for child in value:
                self.add_child(child) # ensure each child is added properly
        else:
            self.add_child(value) # type is not validated, but we assume it's a GFF3Feature object
    
    @property
    def parents(self):
        return self._parents
    
    @parents.setter
    def parents(self, value):
        self._parents = value if isinstance(value, set) \
                        else set(value) if isinstance(value, list) \
                        else set([value]) if isinstance(value, str) \
                        else set()
    
    def add_child(self, childFeature):
        '''
        Adds a child feature to this feature's children list.
        
        Parameters:
            childFeature -- a GFF3Feature object to add as a child of this feature.
        '''
        childFeature.parents.add(self.ID) # ensure the child knows its parent
        self.children.append(childFeature)
        self.__dict__.setdefault(childFeature.ftype, [])
        self.__dict__[childFeature.ftype].append(childFeature)
    
    def find_with_children(self, attribute, foundChildren=None):
        '''
        Finds ("catches") all children contained under this feature which have the specified attribute.
        
        Parameters:
            attribute -- a string indicating an attribute that a child feature should have.
                         For example, to retrieve all "CDS" children at any level under this
                         feature, you would provide "CDS" as the attribute.
            foundChildren -- a list that stores values during recursion, or None for the top-level
                             recursion.
        '''
        if foundChildren is None:
            foundChildren = []
        
        if len(self.children) != 0:
            for child in self.children:
                if hasattr(child, attribute):
                    foundChildren.append(child)
                else:
                    return child.find_with_children(attribute, foundChildren)
        return foundChildren
    
    def find_all_children(self, foundChildren=None):
        '''
        Returns _all_ children contained under this feature.
        '''
        if foundChildren is None:
            foundChildren = []
        
        if len(self.children) != 0:
            for child in self.children:
                foundChildren.append(child)
                child.find_all_children(foundChildren)
        return foundChildren
    
    def length(self, attribute):
        '''
        Sums the length of all .attribute children contained under this feature.
        
        Parameters:
            attribute -- a string indicating the attribute under which children are indexed
        Returns
            length -- the summed length of children
        '''
        if not hasattr(self, attribute):
            return None
        thisLength = 0
        for child in getattr(self, attribute):
            thisLength += (child.end - child.start + 1)
        return thisLength
    
    def format(self, alreadyFound, recursion=None):
        '''
        This method will attempt to render a GFF3-correct format of the
        data this GFF3Feature contains. Several assumptions are made which, if you
        haven't done anything truly weird, will hold true.
        
        alreadyFound should always be empty for the parent-level feature
        
        Parameters:
            alreadyFound -- a set maintained by the calling function which tracks
                            feature .ID values that have already been found/formatted.
                            Necessary to accommodate multiparent features and prevent
                            duplicate outputs.
        '''
        # Recursion management
        if recursion == None:
            recursion = []
        if self.ID in alreadyFound:
            return
        alreadyFound.add(self.ID)
        
        # Depth-first formatting of details
        recursion.append(str(self))
        for child in self.children:
            child.format(alreadyFound, recursion)
        return "\n".join(recursion) + "\n"
    
    def as_sequence(self, fastaObj):
        '''
        Making use of the class objects defined in this repository, obtains
        the corresponding sequence portion identified by this GFF3Feature.
        
        Parameters:
            fastaObj -- a FASTA or Sequence object
        Returns:
            seqObj -- a Sequence object with start and end defined by
                      self.start and self.end (and optionally self.contig)
        '''
        if hasattr(fastaObj, "isFASTA") and fastaObj.isFASTA:
            return fastaObj(self.contig, self.start-1, self.end) # counteract 1-based GFF3 numbering
        elif hasattr(fastaObj, "isSequence") and fastaObj.isSequence:
            return fastaObj[self.start-1:self.end] # counteract as well
        else:
            raise TypeError(f"Cannot obtain '{type(fastaObj)}' type as sequence")
    
    def as_gene_model(self, fastaObj, sequenceType="CDS"):
        '''
        Making use of the class objects defined in this repository, this function
        will order children of CDS or exon type and format this as a Sequence
        object corresponding to the annotation this feature represents.
        '''
        ACCEPTED_TYPES = ["CDS", "exon"]
        if not sequenceType in ACCEPTED_TYPES:
            raise ValueError(f"Feature cannot be ordered by '{sequenceType}'; must be in the list '{ACCEPTED_TYPES}'")
        if not hasattr(self, sequenceType):
            raise ValueError(f"'{self.ID}' lacks any children with '{sequenceType}' type")
        
        # Get the sequence as a string
        startingFrame = None
        sequence = ""
        for subFeature in sorted(getattr(self, sequenceType), key = lambda x: (x.start, x.end)):
            sequence += str(subFeature.as_sequence(fastaObj))
            if startingFrame == None: # capture the first frame for +ve strand features
                startingFrame = subFeature.frame
        if self.strand == "-":
            'Sort is from genomic left->right; first exon of a -ve strand feature is the rightmost'
            startingFrame = subFeature.frame
        
        # Convert to Sequence object
        seqObj = Sequence(self.ID, sequence, startingFrame)
        
        # Reverse complement if necessary
        if self.strand == "-":
            seqObj = seqObj.reverse_complement()
        
        return seqObj
    
    def __str__(self):
        return "\t".join([
            str(getattr(self, x)) if x != "attributes" else
            ";".join([ f"{k}={v}" for k,v in getattr(self, x).items() ])
            for x in GFF3Feature.HEADER_FORMAT
        ])
    
    def __repr__(self):
        reprPairs = []
        attrsToShow = ["ID", "ftype", "contig", "start", "end", "strand", "parents"]
        
        for attr in attrsToShow:
            reprPairs.append("{0}={1}".format(attr, getattr(self, attr)))
        
        return "<{0};{1}>".format(";".join(reprPairs),
                                  f"children=[{', '.join([child.ID for child in self.children])}]")

class UnsortedGff3Error(Exception):
    pass
class DuplicateFeatureError(Exception):
    pass

class GFF3Graph:
    PARENT_INFERENCE = {
        "CDS": "mRNA",
        "exon": "mRNA",
        "mRNA": "gene",
        "lnc_RNA": "gene",
        "Product": "gene" # Product is a special case, but we treat it as a gene parent
        # "gene": None  # Gene is the top-level feature, no parent should be inferred
    }
    SEMICOLON_REGEX = re.compile(r";{2,}")
    
    @staticmethod
    def format_attributes(attributes):
        '''
        Despite how overengineered this may seem, it is necessary to handle cases where multiple redundant
        semi-colons exist, or when e.g., chemical names are embedded which may contain semi colons, not as a 
        delimiter, but as part of the value.
        '''
        # Run an initial cleaning of the string
        attributes = attributes.strip("\r\n\t; ")
        multiSemiColons = sorted(set(GFF3Graph.SEMICOLON_REGEX.findall(attributes)), key=len, reverse=True)
        for multipleSemiColonString in multiSemiColons:
            attributes = attributes.replace(multipleSemiColonString, ";")
        
        # Parse string into pairs of key:value attributes
        splitAttributes = []
        for a in attributes.split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        return { splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2) }
    
    def __init__(self, fileLocation, deduplicate=False):
        self.fileLocation = fileLocation
        self.ftypes = {} # stores ftype:[featureIDs...]
        self.parentFtypes = set() # stores keys in self.ftypes that point to parent-level features
        self.features = {} # stores featureID:GFF3Feature
        self.contigs = set() # stores unique feature .contig values
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        self._nclsMax = None
        
        self.parse_gff3(self.fileLocation, deduplicate)
        self.isGFF3Graph = True # flag for easier type checking
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise ValueError("File location must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"GFF3 file '{value}' is not a file")
        
        self._fileLocation = value
    
    def _format_new_feature_id(self, ftype):
        '''
        Parameters:
            ftype -- a string of a feature type that exists in self.ftypes
        Returns:
            newFeatureID -- a (most likely to be) unique feature ID, although
                            this is not guaranteed for weird GFF3s that are
                            trying to break things as edge cases
        '''
        return f"{ftype}.{len(self.ftypes[ftype]) + 1}"
    
    def _get_unique_feature_id(self, inputID, separator="."):
        ongoingCount = 1
        featureID = inputID
        while featureID in self.features:
            featureID = f"{inputID}{separator}{ongoingCount}"
            if not featureID in self.features:
                break
            ongoingCount += 1
        return featureID
    
    def parse_gff3(self, gff3File, deduplicate=False):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3File -- a string pointing to a GFF3 file location
                        to parse and populate this object with.
            deduplicate -- (OPTIONAL) a boolean indicating whether
                           duplicate feature IDs should raise an error
                           (False) or if we should automatically
                           deduplicate feature IDs (True) to avoid errors
        '''
        # Reset the graph
        self.fileLocation = gff3File
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        # Parse the GFF3 file into a graph structure
        lineCount = 0
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                lineCount += 1
                sl = line.strip("\r\n\t;'\" ").split("\t")
                
                # Skip filler and comment lines
                if line.startswith("#") or len(sl) != 9:
                    continue
                
                # Extract information from this line
                contig, source, ftype, start, end, \
                    score, strand, frame, attributes = sl
                start = int(start)
                end = int(end)
                ftype = GFF3Feature.make_ftype_case_appropriate(ftype)
                attributes = GFF3Graph.format_attributes(attributes)
                
                # Establish or populate tracking containers
                self.ftypes.setdefault(ftype, [])
                self.contigs.add(contig)
                
                # Get the ID attribute
                if not "ID" in attributes:
                    featureID = self._format_new_feature_id(ftype)
                else:
                    featureID = attributes["ID"]
                
                # Get the parent ID(s) attribute
                if not "Parent" in attributes:
                    parentIDs = []
                    self.parentFtypes.add(ftype)
                else:
                    parentIDs = [ x.strip() for x in attributes["Parent"].split(",") ]
                
                # Create a feature object
                feature = GFF3Feature(ID=featureID, ftype=ftype,
                                      start=start, end=end, strand=strand,
                                      source=source, score=score, frame=frame, attributes=attributes,
                                      contig=contig, children=[], parents=parentIDs)
                
                # Index the feature if it doesn't already exist
                if not featureID in self.features:
                    self.add(feature)
                # Specifically handle exons or CDS which are allowed to have duplicated IDs
                elif ftype in ["exon", "CDS"]:
                    feature.ID = self._get_unique_feature_id(featureID)
                    self.add(feature)
                # Handle other duplicated feature types
                else:
                    existingFeature = self.features[featureID]
                    
                    # Handle inferred features
                    if existingFeature.isInferred:
                        # Check that inferred details are correct
                        if existingFeature.ftype != ftype:
                            raise UnsortedGff3Error(f"Feature ID '{featureID}' has a different type '{ftype}' than previously inferred '{existingFeature.ftype}'")
                        if existingFeature.contig != contig:
                            raise UnsortedGff3Error(f"Feature ID '{featureID}' has a different contig '{contig}' than previously inferred '{existingFeature.contig}'")
                        
                        # Update feature details
                        existingFeature.start = start
                        existingFeature.end = end
                        existingFeature.strand = strand
                        existingFeature.source = source
                        existingFeature.score = score
                        existingFeature.frame = frame
                        existingFeature.attributes = attributes
                        existingFeature.parents.update(parentIDs) # add parents to existing set
                        existingFeature.isInferred = False # turn off flag since we found the feature
                        self.add(existingFeature) # re-add to ensure parents are updated correctly
                    
                    # Handle truly duplicated features
                    else:
                        if deduplicate:
                            newFeatureID = self._format_new_feature_id(ftype)
                            feature.ID = newFeatureID # although ID should be "immutable" it is not added into this graph yet so this is acceptable
                            self.add(feature)
                        else:
                            raise DuplicateFeatureError(f"Feature ID '{featureID}' occurs more than once in file '{self.fileLocation}'")
    
    def add(self, feature):
        # Store a new feature within the graph
        if not feature.ID in self.features:
            self.ftypes.setdefault(feature.ftype, []) # necessary if first occurrence of a subfeature preceeds its parent type
            self.ftypes[feature.ftype].append(feature.ID)
            self.features[feature.ID] = feature
        
        # Update graph features with parent-child relationships
        for parentID in feature.parents:
            # Associate the feature with its existing parents
            if parentID in self.features:
                self.features[parentID].add_child(feature)
            # Create a placeholder for the parent if it doesn't exist
            else:
                if feature.ftype in GFF3Graph.PARENT_INFERENCE:
                    parentFeature = GFF3Feature(parentID, GFF3Graph.PARENT_INFERENCE[feature.ftype],
                                                contig=feature.contig, children=[feature],
                                                isInferred=True)
                    self.add(parentFeature)
                else:
                    raise ValueError("Your GFF3 is not sorted in top-down hierarchical order which has caused an error; " +
                                     f"I encountered a {feature.ftype} with ID '{feature.ID}' that has a parent '{parentID}' which has " + 
                                     f"not yet appeared in your GFF3 file. I am unsure what parent type to infer for " +
                                     f"'{feature.ftype}' features, so I cannot continue parsing. Sort your GFF3 file in " +
                                     "conventional top-down hierarchical order before trying again.")
    
    def qc(self, typesToCheck=None):
        '''
        Runs a quality control check on the GFF3Graph object to ensure that all
        features are properly linked.
        
        Prints a warning if any features are found that have no parents or children.
        
        Parameters:
            typesToCheck -- an iterable of strings indicating the feature types to check
                            for dangling features. If None, all feature types are checked.
        '''
        danglingFeatures = {}
        for feature in self:
            if typesToCheck == None or feature.ftype in typesToCheck:
                if len(feature.parents) == 0 and len(feature.children) == 0:
                    danglingFeatures.setdefault(feature.ftype, 0)
                    danglingFeatures[feature.ftype] = 1
        
        if len(danglingFeatures) != 0:
            print(f"WARNING: Parsing '{self.fileLocation}' resulted in dangling features with no parents or children, " +
                  "which is likely due to an unsorted or incorrectly formatted GFF3 file. This may cause issues with " +
                  "BINge's functionality.")
            for ftype, count in danglingFeatures.items():
                print(f"# {count} '{ftype}' feature{'s have' if count > 1 else ' has'} no parents or children")
    
    @property
    def parents(self):
        '''
        Iterates through this GFF3Graph object to yield all features that are
        at the parent level. Beginning with self.parentFtypes, a set which
        only contains ftypes with at least one parent-level feature, we loop
        through all features under those ftype(s) and yield features which
        lack any parents.
        '''
        for ftype in sorted(self.parentFtypes):
            for featureID in self.ftypes[ftype]:
                feature = self[featureID]
                if len(feature.parents) == 0: # we don't trust the GFF3 file to have an ftype ALWAYS be a parent or a child
                    yield feature
    
    def get_feature_parents(self, feature, parents=None):
        '''
        Climbs up any parent(s) to get to the top-level of a given feature. Will return
        the input feature if it is already at the top-level.
        '''
        if parents is None:
            parents = []
        
        if len(feature.parents) == 0:
            parents.append(feature)
        else:
            for parentID in feature.parents:
                return self.get_feature_parents(self[parentID])
        
        return parents
    
    def binge_iterator(self, isMicrobial=False):
        '''
        Makes use of this GFF3Graph class for a specific kind of iteration
        that BINge makes use of. This is being implemented in order to enable
        plug-and-play substitution of the previous GFF3 parsing logic of
        BINge with the newer GFF3Graph approach.
        '''
        # Configuration for microbial or normal GFF3
        if isMicrobial:
            subfeature = ["gene"]
        else:
            subfeature = ["mRNA", "transcript"]
        
        warnedOnce = False
        for subftype in subfeature:
            if subftype in self.ftypes:
                for featureID in self.ftypes[subftype]:
                    feature = self.features[featureID]
                    
                    # Render a CDS coordinates list from the CDS features
                    try:
                        cdsCoords = sorted([ [f.start, f.end, f.frame] for f in feature.CDS ])
                    except:
                        # Handle CDS feature absence via warning
                        if not warnedOnce:
                            print(f"WARNING: '{feature.ID}' lacks CDS features; " + 
                                f"similar warnings for file '{self.fileLocation}' will be suppressed.")
                            warnedOnce = True
                        continue
                    
                    # Render a exon coordinates list from the exon features
                    if isMicrobial:
                        exonCoords = cdsCoords
                    else:
                        # Handle exon feature absence silently
                        "If we get to here, we have a CDS feature, so the exon becomes mostly irrelevant"
                        try:
                            exonCoords = sorted([ [f.start, f.end, f.frame] for f in feature.exon ])
                        except:
                            exonCoords = cdsCoords
                    
                    yield feature.ID, feature.contig, feature.strand, \
                          exonCoords, cdsCoords
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __len__(self):
        return len(self.features)
    
    def __iter__(self):
        return iter(self.features.values())
    
    def __contains__(self, item):
        return item.ID in self.features
    
    def has_key(self, key):
        return key in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3Graph object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.ftypes[key])) for key in self.ftypes.keys()])
        )

class GmapGFF3:
    def __init__(self, fileLocation):
            self.fileLocation = fileLocation
            self.isGmapGFF3 = True # flag for easier type checking
    
    @staticmethod
    def _build_data_to_yield(thisFeature):
        # Parse out variables from mRNA line
        contig, _, _, start, end, \
            _, strand, _, attributes \
            = thisFeature[1]
        
        # Start building a data structure to yield
        dataDict = {"contig": contig, "start": int(start), "end": int(end), "strand": strand}
        
        # Get the GMAP alignment quality attributes
        for attribute in attributes.split(";"):
            key, value = attribute.split("=")
            if key in ["identity", "coverage"]:
                dataDict[key] = float(value)
            elif key == "indels":
                dataDict[key] = int(value)
            elif key == "Name":
                dataDict["Name"] = value
        
        # Get the exons
        dataDict["exons"] = []
        for _sl in thisFeature[2:]:
            contig, _, featureType, start, end, \
                _, strand, _, attributes \
                = _sl
            if featureType == "exon":
                dataDict["exons"].append([int(start), int(end)])
        dataDict["exons"] = sorted(dataDict["exons"])
        
        return dataDict
    
    def __iter__(self):
        '''
        Provides a simple iterator for a GMAP GFF3 file which yields only the
        relevant details for BINge. Assumes the file is sorted as is the case
        with GMAP output.
        '''
        thisAlignment = []
        with open(self.fileLocation, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\t\r\n ").split("\t")
                
                # Skip comment / irrelevant lines
                if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                    continue
                
                # Handle content lines
                else:
                    featureType = sl[2]
                    
                    # Yield a completed gene data dictionary
                    if featureType == "gene" and thisAlignment != []:
                        dataDict = GmapGFF3._build_data_to_yield(thisAlignment)
                        
                        # Reset our feature storage, and yield result now
                        thisAlignment = []
                        yield dataDict
                    
                    # Build an ongoing feature
                    thisAlignment.append(sl)
        
        # Yield the last feature in the GFF3
        if thisAlignment != []:
            dataDict = GmapGFF3._build_data_to_yield(thisAlignment)
        else:
            raise Exception(f"'{self.fileLocation}' does not appear to be a valid GFF3 file!")
        yield dataDict
