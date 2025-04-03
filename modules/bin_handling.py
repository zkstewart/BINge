import os
from .thread_workers import ReturningProcess
from .parsing import load_sequence_length_index
from .gff3_handling import gff3_iterator, iterate_gmap_gff3
from .bins import Bin, BinCollection

# Multithreaded functions and classes
class CollectionSeedProcess(ReturningProcess):
    '''
    Allows for generate_bin_collections() to be run in parallel for multiple genomes.
    
    Parameters:
        gff3File -- a string pointing to a GFF3 file giving annotations for a genome,
                    or None if no pre-seeding is to occur.
        isMicrobial -- a boolean indicating whether the genomes are microbial or not which,
                       in turn, determines whether we will parse mRNA features (False) or
                       gene features (True).
    '''
    def task(self, gff3File, isMicrobial=False):
        binCollection = BinCollection()
        
        # Seed bin collection if GFF3 is available
        if gff3File != None:
            for mrnaID, contig, strand, exon, cds in gff3_iterator(gff3File, isMicrobial):
                # Create a bin for each exon feature
                exonBins = []
                for exonStart, exonEnd, _ in exon:
                    exonBin = Bin(contig, exonStart, exonEnd)
                    exonBin.add(mrnaID)
                    exonBins.append(exonBin)
                
                # Iteratively handle exon bins
                for exonBin in exonBins:
                    binOverlap = find_overlapping_bins(binCollection, exonBin)
                    add_bin_to_collection(binCollection, binOverlap, exonBin)
        
        return binCollection

class GmapBinProcess(ReturningProcess):
    '''
    Bins GMAP alignments into each genome's BinCollection object.
    
    Parameters:
        gmapFile -- a string indicating the location of a GMAP GFF3 for parsing.
        binCollection -- an existing BinCollection object to add GMAP alignments to.
        indexFileName -- a string indicating the location of a pickled dictionary
                         linking sequence IDs (key) to lengths (value).
    '''
    def task(self, gmapFiles, binCollection, indexFileName):
        # Behavioural parameters (static for now, may change later)
        "These statics are for filtering GMAP alignments that are poor quality"
        OKAY_COVERAGE = 96.5
        OKAY_INDEL_PROPORTION = 0.01
        ##
        "These statics prevent a gene mapping to multiple locations when it has a clear best"
        ALLOWED_COV_DIFF = 1
        ALLOWED_IDENT_DIFF = 0.1
        ##
        "This static allows coverage to be lower if the alignment is close to the contig's edge"
        LENIENT_BOUNDARY_LENGTH = 1000
        
        # Load in the lengths index
        seqLenDict = load_sequence_length_index(indexFileName)
        
        # Iterate through GMAP files
        for gmapFile in gmapFiles:
            # Hold onto the GMAP alignments for each gene
            """GMAP is guaranteed to return alignments for each gene together, but not ordered
            by quality. We'll sort them by quality and then process them iteratively"""
            thisBlockID, thisBlockData = None, None
            for dataDict in iterate_gmap_gff3(gmapFile):
                # Handle first iteration
                if thisBlockID == None:
                    thisBlockID = dataDict["Name"]
                    thisBlockData = [dataDict]
                    continue
                
                # Store data if this is from the same sequence
                if dataDict["Name"] == thisBlockID:
                    thisBlockData.append(dataDict)
                    continue
                
                # Otherwise, sort the data dicts by their quality
                "Higher coverage -> higher identity -> lower indels"
                thisBlockData.sort(key = lambda x: (-x["coverage"], -x["identity"], x["indels"]))
                
                # Iterate through data dicts
                bestCov, bestIdent = None, None
                for blockDataDict in thisBlockData:
                    # Get alignment statistics
                    coverage, identity = blockDataDict["coverage"], blockDataDict["identity"]
                    exonLength = sum([
                        end - start + 1
                        for start, end in blockDataDict["exons"]
                    ])
                    indelProportion = blockDataDict["indels"] / exonLength
                    
                    # Apply leniency if the alignment is close to the contig's boundaries
                    """This helps to deal with fragmented genes aligning well to a contig's
                    edge, but the full length model not being binned because its coverage is
                    too low."""
                    beLenient = False
                    if coverage < OKAY_COVERAGE:
                        contigLength = seqLenDict[blockDataDict["contig"]]
                        alignmentStart = min(min(blockDataDict["exons"]))
                        alignmentEnd = max(max(blockDataDict["exons"]))
                        
                        if (alignmentStart < LENIENT_BOUNDARY_LENGTH) \
                            or (alignmentEnd + LENIENT_BOUNDARY_LENGTH) > contigLength:
                            beLenient = True
                    
                    # Skip processing if the alignment sucks
                    isGoodAlignment = (True if beLenient else coverage >= OKAY_COVERAGE) \
                                    and indelProportion <= OKAY_INDEL_PROPORTION
                    if not isGoodAlignment:
                        continue
                    
                    # See if this path should be skipped because another better one was processed
                    if bestCov != None:
                        if (coverage + ALLOWED_COV_DIFF) < bestCov \
                            or (identity + ALLOWED_IDENT_DIFF) < bestIdent:
                            continue
                    else:
                        bestCov, bestIdent = coverage, identity # set if this is the 'best'
                    
                    # Iteratively handle exon features
                    for exonStart, exonEnd in blockDataDict["exons"]:
                        # Create a bin for each exon feature
                        exonBin = Bin(blockDataDict["contig"], exonStart, exonEnd)
                        exonBin.add(blockDataDict["Name"])
                        
                        # See if this overlaps an existing bin
                        binOverlap = find_overlapping_bins(binCollection, exonBin)
                        
                        # Add a new bin, or merge any bins as appropriate
                        add_bin_to_collection(binCollection, binOverlap, exonBin)
                
                # Reset for next iteration
                thisBlockID, thisBlockData = dataDict["Name"], [dataDict]
        
        return binCollection

def add_bin_to_collection(binCollection, binOverlap, newBin):
    '''
    Assistant function to add a new bin on the basis of whatever bins it is
    overlapping. Modifies the binCollection input in place.
    
    Parameters:
        binCollection -- a BinCollection() object from which binOverlap was generated,
                         and the newBin is intended to be added to.
        binOverlap -- a list containing zero or more Bin() objects that the newBin
                      is overlapping.
        newBin -- a Bin() for addition to the binCollection input.
    '''
    # If the newBin does not overlap anything, add it in as-is
    if len(binOverlap) == 0:
        binCollection.add(newBin)
    
    # If it overlaps only one exon, merge it in normally
    elif len(binOverlap) == 1:
        overlappingBin = binOverlap[0]
        
        newBin.merge(overlappingBin)
        binCollection.delete(overlappingBin)
        binCollection.add(newBin)
    
    # It it has multiple overlaps, add its ID in without a true merge
    else:
        for overlappingBin in binOverlap:
            overlappingBin.union(newBin.ids)

def find_overlapping_bins(binCollection, binQuery):
    '''
    Assistant function to help with taking the bins found to overlap a query region
    via the IntervalTree underlying a BinCollection, and then filtering those bins
    to only those that share a substantial proportion of overlap.
    
    Parameters:
        binCollection -- a BinCollection() object to use the .find() method on.
        binQuery -- a Bin() which we want to find overlapping bins from the binCollection
                    from.
    Returns:
        binOverlap -- a list of Bin() objects that overlap the binQuery, and share a
                      substantial proportion of overlap.
    '''
    foundBins = binCollection.find(binQuery.contig, binQuery.start, binQuery.end)
    
    binOverlap = []
    for bin in foundBins:
        if binQuery.is_overlapping(bin):
            binOverlap.append(bin)
    return binOverlap

# Other functions
def generate_bin_collections(genomesDir, threads, isMicrobial):
    '''
    Receives a list of genome FASTA files and generates bin collections each of these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        genomesDir -- a string indicating an existing directory with genome GFF3
                      and/or FASTA files in it.
        threads -- an integer indicating how many threads to run.
        isMicrobial -- a boolean indicating whether the genomes are microbial or not which,
                       in turn, determines whether we will parse mRNA features (False) or
                       gene features (True).
    Returns:
        collectionList -- a list containing BinCollections in numerical order of the genome
                          files.
    '''
    FILE_PREFIX = "genome"
    
    # Locate subdirectory containing files
    if not os.path.isdir(genomesDir):
        raise FileNotFoundError(f"generate_bin_collections failed because '{genomesDir}' doesn't exist or isn't a directory.")
    
    # Locate all genome/GFF3 pairings
    filePairs = []
    for file in os.listdir(genomesDir):
        if file.endswith(".fasta"):
            if not file.startswith(FILE_PREFIX):
                raise ValueError(f"FASTA file '{file}' in '{genomesDir}' does not start with '{FILE_PREFIX}' as expected.")
            
            # Extract file prefix/suffix components
            filePrefix = file.split(".fasta")[0]
            suffixNum = filePrefix.split(FILE_PREFIX)[1]
            if not suffixNum.isdigit():
                raise ValueError(f"FASTA file '{file}' in '{genomesDir}' does not have a number suffix as expected.")
            
            # Add value to pairings list
            filePairs.append([None, suffixNum])
            
            # Add any corresponding GFF3 file if one exists
            gff3File = f"{FILE_PREFIX}{suffixNum}.gff3"
            if os.path.exists(os.path.join(genomesDir, gff3File)):
                filePairs[-1][0] = os.path.join(genomesDir, gff3File)
    
    if not len(filePairs) > 0:
        raise FileNotFoundError(f"generate_bin_collections failed because '{genomesDir}' contains no genome files")
    
    # Sort the list for consistency of ordering
    "If there are gaps in the suffixNum's, ordering is important to keep things paired up"
    filePairs.sort(key = lambda x: int(x[1]))
    suffixes = [ int(x[1]) for x in filePairs ]
    isConsecutive = suffixes == list(range(1, len(suffixes)+1))
    if not isConsecutive:
        raise ValueError(f"generate_bin_collections failed because genome files in '{genomesDir}' are not " + 
                         "consecutively numbered, which is important for later program logic!")
    
    # Start up threads
    collectionList = []
    for i in range(0, len(filePairs), threads): # only process n (threads) at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(filePairs): # parent loop may excess if n > the number of GMAP files
                gff3File, _ = filePairs[i+x]
                
                seedWorkerThread = CollectionSeedProcess(gff3File, isMicrobial) # returns empty BinCollection if no GFF3
                seedWorkerThread.start()
                processing.append(seedWorkerThread)
        
        # Wait on processes to end
        for seedWorkerThread in processing:
            binCollection = seedWorkerThread.get_result()
            seedWorkerThread.join()
            seedWorkerThread.check_errors()
            
            collectionList.append(binCollection)
    
    return collectionList

def populate_bin_collections(genomesDir, collectionList, gmapFiles,
                             threads):
    '''
    Receives a list of BinCollection objects, alongside a list of GMAP GFF3 file
    locations, and uses multiple threads to parse the GMAP files and add them into
    an existing bin, or into a novel bin collection which is returned.
    
    Parameters:
        genomesDir -- a string indicating an existing directory with genome GFF3
                      and/or FASTA files in it.
        collectionList -- a list of BinCollection objects as resulting from
                          generate_bin_collections().
        gmapFiles -- a list of strings pointing to GMAP GFF3 files for population.
        threads -- an integer indicating how many threads to run at a time; this code is
                   parallelised in terms of processing multiple GMAP files at a time,
                   if you have only 1 GMAP file then only 1 thread can be used.
    Returns:
        collectionList -- a list of BinCollections, each BinCollection containing
                          Bins populated with GMAP alignments.
    '''
    # Establish lists for feeding data into threads
    threadData = []
    for suffixNum in range(1, len(collectionList) + 1):
        # Figure out which genome we're looking at
        genomePrefix = f"genome{suffixNum}"
        genomeIndex = int(suffixNum) - 1
        
        # Get the location of the genome length index
        thisGenomeIndex = os.path.join(genomesDir, f"{genomePrefix}.fasta.lengths.pkl")
        
        # Get all GMAP files associated with this genome
        thisGmapFiles = [ gmFile for gmFile in gmapFiles if f"to_{genomePrefix}_" in gmFile]
        
        # Get other data structures for this genome
        thisBinCollection = collectionList[genomeIndex]
        
        # Store for threading
        threadData.append([thisGmapFiles, thisBinCollection, thisGenomeIndex])
    
    # Start up threads
    collectionList = []
    for i in range(0, len(threadData), threads): # only process n (threads) collections at a time
        processing = []
        for x in range(threads): # begin processing n collections
            if i+x < len(threadData): # parent loop may excess if n > the number of GMAP files
                thisGmapFiles, thisBinCollection, thisGenomeIndex = threadData[i+x]
                
                populateWorkerThread = GmapBinProcess(thisGmapFiles, thisBinCollection,
                                                      thisGenomeIndex)
                populateWorkerThread.start()
                processing.append(populateWorkerThread)
        
        # Wait on processes to end
        for populateWorkerThread in processing:
            binCollection = populateWorkerThread.get_result()
            populateWorkerThread.join()
            populateWorkerThread.check_errors()
            
            collectionList.append(binCollection)
        
    return collectionList
