#! python3
# BINge.py
# BIN Genes for Expression analyses

# This script/program aims to provide an alternative to Corset
# for situations where a reference genome is available. However,
# in the reference genome might be for a different subspecies than
# the one you're working with. Hence, it's not good enough to
# just map against the reference genome. You should map against a
# de novo transcriptome, but leverage the genomic information
# to group transcripts into genes. That's what this does.

import os, argparse, sys
from Bio import SeqIO
from intervaltree import IntervalTree

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts import ZS_GFF3IO

# Define classes
class Bin:
    '''
    Make sure to use this bin with 1-based indexing for start and end values,
    in the same way as a GFF3 would.
    '''
    def __init__(self, contig, start, end, representative, binType):
        self.contig = contig
        self.start = start
        self.end = end
        self.representative = representative # should be a feature
        self.binType = binType
        
        self.ids = set() # set of sequence IDs
        self.links = [] # list of other Bin objects this one connects to
        
        self.add_id(representative.ID)
    
    @property
    def contig(self):
        return self._contig
    
    @contig.setter
    def contig(self, value):
        assert isinstance(value, str)
        assert value != ""
        self._contig = value
    
    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self, value):
        assert isinstance(value, int)
        assert value > 0, \
            "A Bin must start at a value >= 1 (it behaves 1-based for indexing)"
        self._start = value
    
    @property
    def end(self):
        return self._end
    
    @end.setter
    def end(self, value):
        assert isinstance(value, int)
        assert value > 0, \
            "A Bin must end at a value >= 1 (it behaves 1-based for indexing)"
        self._end = value
    
    @property
    def representative(self):
        return self._representative
    
    @representative.setter
    def representative(self, value):
        self._representative = value
    
    @property
    def binType(self):
        return self._binType
    
    @binType.setter
    def binType(self, value):
        self._binType = value
    
    def add_id(self, idValue):
        self.ids.add(idValue)
    
    def add_link(self, binValue):
        self.links.append(binValue)
    
    def __repr__(self):
        return (f"<Bin object;binType='{self.binType}';contig='{self.contig}';start={self.start};" +
                f"end={self.end};representative='{self.representative}';num_ids={self.ids}>"
        )

class BinCollection:
    '''
    Encapsulates an IntervalTree allowing easy addition and finding of Bin objects.
    '''
    def __init__(self):
        self.bins = IntervalTree()
    
    def add(self, bin):
        # assert isinstance(bin, Bin) # doesn't work
        self.bins[bin.start:bin.end+1] = bin
    
    def find(self, contig, start, end):
        '''
        Start and end should be provided as 1-based and inclusive values.
        For example, searching for 100,100 will find overlaps at that exact
        position.
        '''
        bins = [ b.data for b in self.bins[start:end+1] if b.data.contig == contig ]
        return bins
    
    def __repr__(self):
        return "<BinCollection object;num_bins={0}>".format(
            len(self.bins)
        )

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.transcriptomeFile):
        print(f'I am unable to locate the transcriptome FASTA file ({args.transcriptomeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for genomeFile in args.genomeFiles:
        if not os.path.isfile(genomeFile):
            print(f'I am unable to locate the genome FASTA file ({genomeFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for annotationFile in args.annotationFiles:
        if not os.path.isfile(annotationFile):
            print(f'I am unable to locate the genome GFF3 file ({annotationFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for gmapFile in args.gmapFiles:
        if not os.path.isfile(gmapFile):
            print(f'I am unable to locate the GMAP GFF3 file ({gmapFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate the input files are logically sound
    if len(args.genomeFiles) != len(args.gmapFiles):
        print("Your genome and GMAP files are incompatible!")
        print(f"I'm seeing {len(args.genomeFiles)} genome files and {len(args.gmapFiles)} GMAP files")
        print("These numbers should be the same. You need to fix this up and try again.")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def _create_feature_from_sl(sl):
    # Extract details from sl
    contig, source, featureType, start, end, \
        score, strand, frame, attributes \
        = sl
    
    # Build feature with main details
    feature = ZS_GFF3IO.Feature()
    feature.add_attributes({
        "contig": contig, "source": source, "type": featureType,
        "start": int(start), "end": int(end), "coords": [int(start), int(end)],
        "score": score, "strand": strand, "frame": frame
    })
    
    # Add attributes
    splitAttributes = []
    for a in attributes.split("="):
        if ";" in a:
            splitAttributes += a.rsplit(";", maxsplit=1)
        else:
            splitAttributes.append(a)
    attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
    feature.add_attributes(attributesDict)
    
    return feature

def iterate_through_gmap(gmapFile):
    '''
    Provides a simple iterator for a GMAP GFF3 file which yields GFF3
    features.
    
    Parameters:
        gmapFile --
        recordsDict -- 
        gff3Obj -- 
        binCollection -- 
    '''
    thisFeature = []
    with open(gmapFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\t\r\n ").split("\t")
            
            # Skip comment / irrelevant lines
            if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                continue
            
            # Handle content lines
            else:
                featureType = sl[2]
                
                # Yield a completed gene feature
                if featureType == "gene" and thisFeature != []:
                    # Create the base gene feature
                    feature = _create_feature_from_sl(thisFeature[0])
                    
                    # Add subfeatures to base
                    for _sl in thisFeature:
                        _feature = _create_feature_from_sl(_sl)
                        feature.add_child(_feature)
                    
                    # Reset our feature storage, and yield result now
                    thisFeature = []
                    yield feature
                
                # Build an ongoing feature
                thisFeature.append(sl)

def calculate_overlap_percentages(geneFeature, feature):
    # Calculate overlap for CDS
    featureCoords = [ pos for coords in ZS_GFF3IO.GFF3._get_feature_coords(feature, "CDS")[0] for pos in coords ]
    featureStart, featureEnd = min(featureCoords), max(featureCoords)
    
    geneCoords = []
    for mrnaFeature in geneFeature.mRNA:
        for CDSFeature in mrnaFeature.CDS:
            geneCoords += [CDSFeature.start, CDSFeature.end]
    geneStart, geneEnd = min(geneCoords), max(geneCoords)
    
    cds_overlapLength = min(geneEnd, featureEnd) - max(geneStart, featureStart)
    cds_geneOvlPct = cds_overlapLength / (geneEnd - geneStart + 1)
    cds_featureOvlPct = cds_overlapLength / (featureEnd - featureStart + 1)
    
    # Calculate overlap for overall
    overlapLength = min(geneFeature.end, feature.end) - max(geneFeature.start, feature.start)
    geneOvlPct = overlapLength / (geneFeature.end - geneFeature.start + 1)
    featureOvlPct = overlapLength / (feature.end - feature.start + 1)
    
    return max(geneOvlPct, cds_geneOvlPct), max(featureOvlPct, cds_featureOvlPct) # be optimistic

def bin_by_gmap(gmapFile, recordsDict, gff3Obj, binCollection):
    '''
    Parameters:
        gmapFile --
        recordsDict -- 
        gff3Obj -- 
        binCollection -- 
    '''
    # Behavioural parameters (static for now, may change later)
    VERY_GOOD_COVERAGE = 99
    VERY_GOOD_IDENTITY = 94
    ##
    GOOD_COVERAGE = 95
    GOOD_IDENTITY = 90
    ##
    OKAY_COVERAGE = 88
    OKAY_IDENTITY = 87
    ##
    ALLOWED_COV_DIFF = 3
    ALLOWED_IDENT_DIFF = 2
    ##
    GOOD_OVL_PCT = 90 # needs investigation
    
    # Parse GMAP file and begin binning
    cov = []
    ident = []
    
    pathDict = {} # holds onto sequences we've built a path from already
    for feature in iterate_through_gmap(gmapFile):
        # Get alignment statistics
        coverage, identity = float(feature.mRNA[0].coverage), float(feature.mRNA[0].identity)
        cov.append(coverage)
        ident.append(identity)
        
        # Skip processing if the alignment sucks
        if not (coverage >= OKAY_COVERAGE and identity >= OKAY_IDENTITY):
            continue
        isGoodAlignment = coverage >= GOOD_COVERAGE and identity >= GOOD_IDENTITY
        
        # See if this path should be skipped because another better one was processed
        baseID = feature.ID.rsplit(".", maxsplit=1)[0]
        if baseID in pathDict:
            prevCov, prevIdent = pathDict[baseID]
            
            # Skip if the first path was the best
            if (coverage + ALLOWED_COV_DIFF) < prevCov \
                or (identity + ALLOWED_IDENT_DIFF) < prevIdent:
                continue
        
        # See if this overlaps an existing gene and/or bin
        geneOverlap = gff3Obj.ncls_finder(feature.start, feature.end, "contig", feature.contig)
        binOverlap = binCollection.find(feature.contig, feature.start, feature.end)
        stop
        
        # GENE CASE
        if len(geneOverlap) == 1:
            geneFeature = geneOverlap[0]
            
            # Calculate how much the projected feature overlaps the gene feature
            geneOvlPct, featureOvlPct = calculate_overlap_percentages(geneFeature, feature)
            isGoodOvl = geneOvlPct >= GOOD_OVL_PCT and featureOvlPct >= GOOD_OVL_PCT
            
            # Subcase: create new bin
            if len(binOverlap) == 0:
                # Subcase: new bin based on a gene
                if coverage >= GOOD_COVERAGE and identity >= GOOD_IDENTITY and isGoodOvl:
                    
                    geneBin = Bin(geneFeature.contig, geneFeature.start,
                                    geneFeature.end, geneFeature, "gene")
                    binCollection.add(geneBin)
                    
                    # Add our transcript to it
                    geneBin.add_id(baseID)
                
                # Subcase: novel bin based on this feature
                elif coverage >= GOOD_COVERAGE and identity >= GOOD_IDENTITY:
                    novelBin = Bin(feature.contig, feature.start,
                                 feature.end, feature, "novel")
                    binCollection.add(novelBin)
                
                # Subcase: 
            
            # Subcase: add to existing bin
            else:
                pass
        
        # NOVEL CASE
        if len(geneOverlap) == 0:
            # Subcase: new bin
            
            # Subcase: add to bin
            
            # Subcase: expand bin
            pass
        
        # EASY CASE
        if len(geneOverlap) == 1 and len(binOverlap) <= 1:
            
            # EASY - VERY GOOD or GOOD
            if coverage >= GOOD_COVERAGE and identity >= GOOD_IDENTITY:
                pass
                # Subcase: new bin based on a gene
                
                
                # Subcase: add to existing gene bin
                # else:
                #     geneBin = binOverlap[0]
                #     geneBin.add_id(baseID)
                    
                #     assert geneBin.binType == "gene", \
                #         "testing assert, how can EASY VERY GOOD get here?"
            
            # EASY - OKAY AND WORSE
            else:
                
                    
                # Subcase: add to existing novel bin
                pass
        
        # EASY OVERLAP CASE
        if len(geneOverlap) == 1 and len(binOverlap) > 1:
            
            # Subcase: add to existing gene bin
            if coverage >= GOOD_COVERAGE and identity >= GOOD_IDENTITY:
                pass
            
            # Subcase: add to existing novel bin
            
            # Subcase: new bin based on
            pass
        
        # CHIMERA CASE
        if len(geneOverlap) >= 2:
            pass
        
        # LOOSE FIT 1
        
        
        # LOOSE FIT 2
        
        # ABSENT CASE
        
        
        # UNHANDLED BASKET
        pass
        
        # Add to path dict (if necessary)
        pathDict.setdefault(baseID, [coverage, identity])

## Main
def main():
    # User input
    usage = """%(prog)s is ... WIP.
    
    For each genome file given as input, you should provide a matching annotation
    GFF3, and a matching GMAP GFF3 alignment file. These values should be ordered
    equivalently.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="transcriptomeFile",
                   required=True,
                   help="Input transcriptome FASTA file (mRNA or CDS).")
    p.add_argument("-ge", dest="genomeFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more genome FASTA files")
    p.add_argument("-ga", dest="annotationFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more genome annotation (GFF3) files")
    p.add_argument("-gm", dest="gmapFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more GMAP (GFF3) files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each GFF3 into a slim structure with NCLS indexing
    gff3s = []
    for annotFile in args.annotationFiles:
        gff3Obj = ZS_GFF3IO.GFF3(annotFile, strict_parse=True, slim_index=True)
        gff3Obj.create_ncls_index(typeToIndex="gene")
        gff3s.append(gff3Obj)
    
    # For each genome and GMAP file, add sequences to bins
    binCollection = BinCollection()
    for genomeFile, gmapFile, gff3Obj in zip(args.genomeFiles, args.gmapFiles, gff3s):
        # Parse the genomeFile as dict for unsorted lookups
        recordsDict = SeqIO.to_dict(SeqIO.parse(genomeFile, "fasta"))
        
        # Loop through gmap alignments and bin sequences
        stop
        bin_by_gmap(gmapFile, recordsDict, gff3Obj, binCollection)
    
    ## TBD: The whole entirety of the code
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
