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

import os, argparse, sys, queue, pickle
import networkx as nx
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Various_scripts.Function_packages import ZS_ClustIO, ZS_BlastIO

from modules.bins import Bin, BinCollection
from modules.thread_workers import GmapBinThread, BinSplitWorkerThread, CollectionWorkerThread
from modules.gff3_handling import GFF3
from modules.fasta_handling import ZS_SeqIO, FastaCollection

# Define functions
def validate_args(args):
    # Validate input file locations
    for fastaFile in args.fastaFiles:
        if not os.path.isfile(fastaFile):
            print(f'I am unable to locate the transcript FASTA file ({fastaFile})')
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
    if len(args.annotationFiles) != 0:
        if len(args.annotationFiles) != len(args.gmapFiles):
            print("Your genome annotation and GMAP files are incompatible!")
            print(f"I'm seeing {len(args.annotationFiles)} annotation files and {len(args.gmapFiles)} GMAP files")
            print("These numbers should be the same. You need to fix this up and try again.")
            quit()
    else:
        print("Note that you are choosing to run BINge without an input annotation")
        print("That's okay if there isn't one available, but if there IS one available, " + 
              "I'd recommend that you use it.")
        print("Program will continue operation as usual.")
    
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
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def generate_bin_collections(annotationFiles):
    '''
    Receives a list of genome annotations in GFF3 format and parses these
    into BinCollection structures which are separately stored in the returned list.
    
    Parameters:
        annotationFiles -- a list of strings pointing to GFF3 genome annotation files.
    Returns:
        collectionList -- a list containing BinCollections in the same order as the
                          input annotation files.
    '''
    # Parse each GFF3 into a bin collection structure
    collectionList = [] # by using a list, we keep genome bins separate to multi-thread later
    
    for annotFile in annotationFiles:
        binCollection = BinCollection()
        
        # Parse first as a GFF3 object
        gff3Obj = GFF3(annotFile, strict_parse=False)
        for geneFeature in gff3Obj.types["gene"]:
            
            # Create a bin for this feature
            featureBin = Bin(geneFeature.contig, geneFeature.start, geneFeature.end)
            try:
                featureBin.add(geneFeature.ID, Bin.format_exons_from_gff3_feature(geneFeature))
            except:
                "This exception occurs if a gene feature has non-mRNA children e.g., ncRNAs"
                continue
            
            # See if this overlaps an existing bin
            binOverlap = binCollection.find(geneFeature.contig, geneFeature.start, geneFeature.end)
            
            # If not, add the new bin
            if len(binOverlap) == 0:
                binCollection.add(featureBin)
            
            # Otherwise...
            else:
                # ... merge the bins together
                for overlappingBin in binOverlap:
                    featureBin.merge(overlappingBin)
                
                # ... delete the overlapping bins
                for overlappingBin in binOverlap:
                    binCollection.delete(overlappingBin)
                
                # ... and add the new bin
                binCollection.add(featureBin)
        
        # Store the bin collection in our list for multi-threading later
        collectionList.append(binCollection)
    
    return collectionList

def populate_bin_collections(collectionList, gmapFiles, threads):
    '''
    Receives a list of BinCollection objects, alongside a matching list
    of GMAP GFF3 file locations, and uses multiple threads to parse the GMAP
    files and add them into an existing bin, or into a novel bin collection
    which is returned.
    
    Parameters:
        collectionList -- a list of BinCollection objects as resulting from
                          generate_bin_collections().
        gmapFiles -- a list of strings pointing to GMAP GFF3 files; they should
                     be ordered the same as the files which were used as input to
                     generate_bin_collections().
        threads -- an integer indicating how many threads to run; this code is
                   parallelised in terms of processing multiple GMAP files at a time,
                   if you have only 1 GMAP file then only 1 thread will be used.
    Modifies:
        collectionList -- the BinCollection objects in this list will have alignment IDs
                          added to the contained Bin objects.
    Returns:
        novelBinCollection -- a BinCollection containing Bins which did not overlap existing
                              Bins part of the input collectionList.
    '''
    novelBinCollection = BinCollection()
    multiOverlaps = []
    for i in range(0, len(gmapFiles), threads): # only process n (threads) files at a time
        processing = []
        for x in range(threads): # begin processing n files
            if i+x < len(gmapFiles): # parent loop may excess if n > the number of GMAP files
                gmapFile = gmapFiles[i+x]
                binCollection = collectionList[i+x]
                
                gmapWorkerThread = GmapBinThread(gmapFile, binCollection)
                processing.append(gmapWorkerThread)
                gmapWorkerThread.start()
        
        # Gather results
        for gmapWorkerThread in processing:
            # Wait for thread to end
            gmapWorkerThread.join()
            
            # Grab the new outputs from this thread
            "Each thread modifies a BinCollection part of collectionList directly"
            threadNovelBinCollection = gmapWorkerThread.novelBinCollection
            multiOverlaps.append(gmapWorkerThread.multiOverlaps)
            
            # Merge them
            novelBinCollection.merge(threadNovelBinCollection)
    return novelBinCollection, multiOverlaps

def iterative_bin_self_linking(binCollection, convergenceIters):
    '''
    Links the bins in a BinCollection to themselves iteratively until
    it converges or convergenceIters is reached.
    
    Parameters:
        binCollection -- a BinCollection object
        convergenceIters -- an integer providing a maximum limit for
                            the amount of iterations that may occur.
    '''
    prevCollectionCount = len(binCollection.bins)
    for _ in range(convergenceIters):
        binCollection = binCollection.link_bins()
        if len(binCollection.bins) == prevCollectionCount:
            break
        prevCollectionCount = len(binCollection.bins)
    return binCollection

def multithread_bin_splitter(binCollection, threads):
    '''
    Using a queue system, worker threads will grab bins from the binCollection parameter
    process with the BinSplitter class. They will put their results into the output
    worker which will help to format a new BinCollection object as output.
    
    We want to split bins because each bin is only guaranteed to contain
    sequences that generally align well over a predicted genomic region. There's no
    guarantee that they are isoforms of the same gene, or even share any exonic content.
    This splitting step is very lax and designed specifically for separating out nested
    genes within the introns of what might be considered the "main gene" the
    bin represents.
    
    Parameters:
        binCollection -- a BinCollection containing Bins which have aligned against
                         a reference genome, and hence have informative .exons values.
        threads -- an integer indicating how many threads to run.
        clusterType -- a string to put in an output column indicating what source of
                       evidence led to this cluster's formation; usually, this will
                       have the value of "binned" to indicate that it was binned via
                       GMAP alignment, or "unbinned" to indicate that it has been
                       binned solely via CD-HIT clustering.
    Returns:
        newBinCollection -- a new BinCollection object to replace the original one.
    '''
    # Set up queueing system for multi-threading
    workerQueue = queue.Queue(maxsize=int(threads * 10)) # allow queue to scale with threads
    outputQueue = queue.Queue(maxsize=int(threads * 20))
    
    # Start up threads for clustering of bins
    workers = []
    for _ in range(threads):
        worker = BinSplitWorkerThread(workerQueue, outputQueue)
        worker.daemon = True
        worker.start()
        workers.append(worker)
    
    outputWorker = CollectionWorkerThread(outputQueue)
    outputWorker.daemon = True
    outputWorker.start()
    
    # Put bins in queue for worker threads
    for interval in binCollection:
        bin = interval.data
        workerQueue.put(bin)
    
    # Close up shop on the threading structures
    for i in range(threads):
        workerQueue.put(None) # this is a marker for the worker threads to stop
    for worker in workers:
        worker.join() # need to call .join() on the workers themselves, not the queue!
    
    outputQueue.put([None, None]) # marker for the output thead to stop; needs 2 Nones!
    outputWorker.join()
    
    return outputWorker.binCollection

def bin_chimera_splitter(binCollection, multiOverlap):
    '''
    Receives a BinCollection that has potentially been created through multiple
    genome, multiple subspecies transcriptomics. Despite this confusion, it will
    attempt to split bins within genomes that are likely to be chimeric.
    
    It does this using 'multiOverlap', a list of GMAP alignment Features which
    overlap more than one Bin. It is assumed that these Bins have been previously
    seeded from the genome's annotation, and hence this function will help us to
    identify situations where the genome annotation has flaws.
    
    The method of performing this is similar to that seen in bin_self_linker()
    where it's leveraging the graph-based approach to model bins which should
    merge. Merging decisions are made based on comparing the 'weight' of an edge
    (i.e., how many GMAP alignments support the bin's merging) to the number of
    IDs which do not support the merging (i.e., how many GMAP alignments only
    overlapped a single bin). The vote threshold decides whether we should merge
    or not.
    
    Note: this function is not guaranteed to produce a stable output i.e., if
    you feed its result back in, you may get different results. It should
    eventually converge on a stable solution, but there might theoretically
    be an edge case where that never happens.
    
    Parameters:
        binCollection -- a BinCollection object containing an assortment
                         of bins from different genomes and/or genes that
                         are (at a sequence level) indistinguishable and
                         hence confounding for DGE.
        multiOverlap -- a list containing GFF3 Features which were found to
                        overlap more than one Bin in the binCollection.
    Returns:
        linkedBinCollection -- a BinCollection object where bins have been
                               merged where deemed appropriate.
    '''
    VOTE_THRESHOLD = 0.5
    
    # Format bins into a dictionary
    binDict = {}
    numBins = 0
    for bin in binCollection:
        binDict[bin.data.sha256()] = bin.data
        numBins += 1
    assert len(binDict) == numBins, \
        "Hash collision occurred! I can't handle this, and you should buy a lottery ticket..."
    
    # Model links between bins as a graph structure
    binGraph = nx.Graph()
    binGraph.add_nodes_from(binDict.keys())
    
    # Process bins that are multi-overlapped
    for overlappingFeature in multiOverlap:
        # Find the bins which are overlapped
        binOverlap = binCollection.find(overlappingFeature.contig, overlappingFeature.start, overlappingFeature.end)
        assert len(binOverlap) > 1, "Found a non-multi-overlapper somehow?"
        binOverlap = [ bin.sha256() for bin in binOverlap ]
        
        # Add edges between overlapped bins
        for i in range(0, len(binOverlap)-1):
            for x in range(i+1, len(binOverlap)):
                baseID = overlappingFeature.ID.rsplit(".", maxsplit=1)[0]
                exons = Bin.format_exons_from_gff3_feature(overlappingFeature)
                
                # Add new edge
                if not binGraph.has_edge(binOverlap[i], binOverlap[x]):
                    binGraph.add_edge(binOverlap[i], binOverlap[x],
                                    weight=1, ids=[baseID],
                                    exons={baseID: exons})
                
                # Increase weight of existing edge
                else:
                    binGraph[binOverlap[i]][binOverlap[x]]["weight"] += 1
                    binGraph[binOverlap[i]][binOverlap[x]]["ids"].append(baseID)
                    binGraph[binOverlap[i]][binOverlap[x]]["exons"][baseID] = exons
    
    # Create a new network with edges set if they exceed the unlinked number of IDs
    finalGraph = nx.Graph()
    finalGraph.add_nodes_from(binDict.keys())
    
    for connectedBins in nx.connected_components(binGraph):
        connectedBins = list(connectedBins)
        
        for i in range(0, len(connectedBins)-1):
            for x in range(i+1, len(connectedBins)):
                binHash1 = connectedBins[i]
                binHash2 = connectedBins[x]
                edgeWeight = binGraph[binHash1][binHash2]["weight"]
                distanceWeight = sum([
                    len(binDict[binHash1].ids), len(binDict[binHash2].ids)
                ])
                
                mergeVote = edgeWeight / (edgeWeight + distanceWeight)
                if mergeVote >= VOTE_THRESHOLD:
                    finalGraph.add_edge(binHash1, binHash2,
                                        ids = binGraph[binHash1][binHash2]["ids"],
                                        exons = binGraph[binHash1][binHash2]["exons"])
    
    # Merge bins on the basis of link identification
    linkedBinCollection = BinCollection()
    for connectedBins in nx.connected_components(finalGraph):
        connectedBins = list(connectedBins)
        
        newBin = binDict[connectedBins[0]]
        binHash1 = connectedBins[0]
        for binHash2 in connectedBins[1:]:
            newBin2 = Bin(newBin.contig, binDict[binHash2].start, binDict[binHash2].end)
            newBin2.ids = set(finalGraph[binHash1][binHash2]["ids"]).union(binDict[binHash2].ids)
            
            newBin2.exons = finalGraph[binHash1][binHash2]["exons"]
            newBin2.exons.update(binDict[binHash2].exons)
            
            newBin.merge(newBin2)
        linkedBinCollection.add(newBin)
    
    return linkedBinCollection

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

def get_unbinned_sequence_ids(binCollectionList, transcriptRecords):
    '''
    Compares one or more BinCollection objects against the transcript sequences
    to see if any sequences indicated in transcriptRecords do not exist in
    any Bins.
    
    Parameters:
        binCollectionList -- a list containing BinCollection's
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences
    Returns:
        unbinnedIDs -- a set containing string value for sequence IDs that were not
                       binned by BINge's main clustering process.
    '''
    binnedIDs = []
    for binCollection in binCollectionList:
        for interval in binCollection:
            bin = interval.data
            binnedIDs.extend(bin.ids)
    binnedIDs = set(binnedIDs)
    
    unbinnedIDs = set()
    for record in transcriptRecords:
        if record.name not in binnedIDs:
            unbinnedIDs.add(record.name)
    
    return unbinnedIDs

def cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args):
    '''
    Runs CD-HIT on the unbinned sequences in order to assign them to a cluster.
    It's not ideal, but the alternative is to exclude these sequences which may
    not be in line with the user's aims. By noting the source of clustering in
    the output file, they can decide if they'd like to exclude these or not.
    
    Note: the default CD-HIT parameter values have been derived from some limited
    testing which suggests that these values may work optimally to cluster sequences
    into "genes" or the closest thing we have to them without proper genomic evidence.
    
    Parameters:
        unbinnedIDs -- a set containg sequence IDs that exist in transcriptRecords
                       which should be clustered with an external algorithm.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences.
        args -- an argparse ArgumentParser object with attributes as set by BINge's
                main argument parsing process.
    '''
    # Generate a temporary FASTA file containing unbinned transcripts
    tmpFileName = write_unbinned_fasta(unbinnedIDs, transcriptRecords)
    
    # Cluster the unbinned transcripts depending on BINge parameters
    if args.unbinnedClusterer in ["mmseqs-cascade", "mmseqs-linclust"]:
        resultClusters = mmseqs_clustering(tmpFileName, args.unbinnedClusterer, 
                                           args.mmseqsDir, args.tmpDir,
                                           args.threads, args.evalue, args.identity,
                                           args.coverage, args.mode, args.sensitivity,
                                           args.steps)
    else:
        resultClusters = cdhit_clustering(tmpFileName, args.cdhitDir, args.threads,
                                          args.cdhitMem, args.identity,
                                          args.cdhitShortCov, args.cdhitLongCov)
    
    # Clean up temporary file
    os.unlink(tmpFileName)
    
    # Return cluster dictionary results
    return resultClusters

def mmseqs_clustering(fastaFile, algorithm, mmseqsDir, tmpDir, threads, evalue, identity,
                      coverage, mode, sensitivity, steps):
    '''
    Parameters:
        fastaFile -- a FASTA file containing nucleotide sequences for clustering.
        algorithm -- a string in the list ["mmseqs-cascade", "mmseqs-linclust"] indicating
                     which clustering algorithm should be used.
        mmseqsDir -- a string indicating the location where the mmseqs executable is found.
        tmpDir -- a string location for where MMseqs2 should keep temp files.
        threads -- a positive integer for how many threads to use when running MMseqs2
                   clustering.
        evalue -- a positive float with a minimum of 0.0 controlling the E-value threshold
                  for clustering.
        identity -- a positive float in the range 0.0 -> 1.0 controlling the sequence identity
                    threshold for clustering.
        coverage -- a positive float in the range 0.0 -> 1.0 controlling the amount of aligned
                    residues in both shorter and longer sequences.
        mode -- a string in the list ["set-cover", "connected-component", "greedy"],
                corresponding to modes 0, 1, and 2,3 of MMseqs2.
    '''
    # Generate the MMseqs2 sequence database
    mmDB = ZS_BlastIO.MM_DB(fastaFile, mmseqsDir, tmpDir, threads)
    mmDB.generate()
    mmDB.index()
    
    # Cluster the unbinned transcripts
    if algorithm == "mmseqs-cascade":
        clusterer = ZS_ClustIO.MM_Cascade(
            mmDB, evalue, identity, coverage,
            mode, threads, tmpDir,
            sensitivity, steps
        )
    else:
        clusterer = ZS_ClustIO.MM_Linclust(
            mmDB, evalue, identity, coverage,
            mode, threads, tmpDir
        )
    clusterer.cluster()
    
    # Generate the tabular output
    tmpFileName = "tmp_BINge_mms2clusttable_{0}.tsv".format(
        ZS_SeqIO.Conversion.get_hash_for_input_sequences(fastaFile)
    )
    clusterer.tabulate(tmpFileName)
    
    # Parse it into a form that BINge can use
    resultClusters = clusterer.parse_tsv(tmpFileName)
    
    # Clean up temporary files
    os.unlink(tmpFileName) # clean up tabular output
    mmDB.clean_all() # clean up sequence database generation and indexing
    clusterer.clean_all() # clean up clustering outputs
    
    # Return cluster dictionary results
    return resultClusters

def cdhit_clustering(fastaFile, cdhitDir, threads, mem,
                     identity, shorterCovPct, longerCovPct, molecule="nucleotide"):
    '''
    Runs CD-HIT on the unbinned sequences in order to assign them to a cluster.
    It's not ideal, but the alternative is to exclude these sequences which may
    not be in line with the user's aims. By noting the source of clustering in
    the output file, they can decide if they'd like to exclude these or not.
    
    Note: the default CD-HIT parameter values have been derived from some limited
    testing which suggests that these values may work optimally to cluster sequences
    into "genes" or the closest thing we have to them without proper genomic evidence.
    
    Parameters:
        fastaFile -- a FASTA file containing nucleotide sequences for clustering.
        cdhitDir -- a string indicating the location where cd-hit-est is found.
        threads -- an integer value indicating how many threads to run CD-HIT with.
        mem -- an integer value indicating how many megabytes of memory to run CD-HIT with.
        identity -- a float value indicating what identity value to run CD-HIT with.
        shorterCovPct -- a float value setting the -aS parameter of CD-HIT.
        longerCovPct -- a float value setting the -aL parameter of CD-HIT.
        molecule -- a string of "nucleotide" or "protein" indicating what molecule
                    the FASTA file sequences are.
    '''
    assert molecule == "nucleotide" or "protein", \
        f"{molecule} must be 'nucleotide' or 'protein'!"
    
    # Cluster the unbinned transcripts
    clusterer = ZS_ClustIO.CDHIT(fastaFile, molecule, cdhitDir)
    clusterer.identity = identity
    clusterer.set_shorter_cov_pct(shorterCovPct)
    clusterer.set_longer_cov_pct(longerCovPct)
    clusterer.set_local()
    clusterer.threads = threads
    clusterer.mem = mem
    clusterer.get_cdhit_results(returnFASTA=False, returnClusters=True)
    
    # Return cluster dictionary results
    return clusterer.resultClusters

def write_unbinned_fasta(unbinnedIDs, transcriptRecords):
    '''
    A helper function which creates a temporary FASTA file containing all unbinned
    sequences. This file will be used for clustering, after which it can be deleted.
    
    Parameters:
        unbinnedIDs -- a 
    '''
    # Generate a temporary FASTA file containing unbinned transcripts
    tmpFileName = "tmp_BINge_unbinned_{0}.fasta".format(
        ZS_SeqIO.Conversion.get_hash_for_input_sequences(str(transcriptRecords))
    )
    with open(tmpFileName, "w") as fileOut:
        for seqID in unbinnedIDs:
            record = transcriptRecords[seqID]
            fileOut.write(f">{record.name}\n{str(record)}\n")
    
    return tmpFileName

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
    
    "These hashes are the only ones which behaviourally influence the pre-external clustering"
    HASHING_PARAMS = ["fastaFiles", "annotationFiles", "gmapFiles", "convergenceIters"]
    
    strForHash = ""
    for param in HASHING_PARAMS:
        strForHash += str(args.__dict__[param])
    paramHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return paramHash

## Main
def main():
    showHiddenArgs = '--help-long' in sys.argv
    
    # User input
    usageShort = """Quick notes for the use of %(prog)s: 1) For each annotation GFF3 (-ga)
    you should provide a matching GMAP GFF3 alignment (-gm). 2) One or more inputs can be
    provided with -i, which are looked at internally as a single file; all sequences in the
    -ga files must be provided with -i (i.e., sequences with the same ID as the GFF3 gene IDs
    must be given) as well as any sequences used for GMAP alignment. 3) Unbinned sequences will be
    cascade clustered using MMseqs2 by default which is recommended; you can choose Linclust
    or CD-HIT if desired. 4) Some parameters are hidden in this short help format since their
    defaults are adequate; specify --help-long to see information for those options.
    """
    
    usageLong = """%(prog)s (BIN Genes for Expression analyses) is a program which bins
    de novo-assembled transcripts together on the basis of reference genome alignments.
    This might be necessary when working with multiple subspecies that are expected to
    diverge only slightly from the reference organism. By binning like this, each subspecies
    can have its gene counts compared fairly during DGE, and some of the pitfalls of
    other approaches e.g., CD-HIT are avoided.
    ###
    Prior to running BINge, you need to run GMAP alignment of your de novo transcriptome
    against one of more reference genomes. Running GMAP with parameters '-f 2 -n 6' is 
    suggested (-f 2 is required to have GFF3 format output!).
    ###
    Note 1: For each annotation GFF3 (-ga) you should provide a matching GMAP GFF3 alignment
    file (-gm). These values should be ordered equivalently.
    ###
    Note 2: You may want to provide multiple inputs with -i, one for your transcripts
    and another for the reference sequences from each of your GFF3(s). All sequences
    indicated in your -ga and -gm files need to be locateable in the files given to -i.
    ###
    Note 3: Sequences which do not align against the genome are considered to be "unbinned".
    These will be clustered with MMseqs2 cascaded clustering by default, which is the recommended
    choice. You can use Linclust (also okay) or CD-HIT (potentially very slow) if wanted.
    ###
    Note 4: You're seeing the --help-long format of this message, which means you may want to
    configure the way clustering of unbinned sequences works. Behavioural parameters of
    the algorithms can be tuned here, but the defaults are expected to work most of the time.
    The main exception is CD-HIT's memory utilisation, which probably should be set depending
    on what you have available. tldr; change these if you know what you're doing.
    """
    
    p = argparse.ArgumentParser(description=usageLong if showHiddenArgs else usageShort)
    # Required
    p.add_argument("-i", dest="fastaFiles",
                   required=True,
                   nargs="+",
                   help="Input transcriptome FASTA file(s) (mRNA or CDS).")
    p.add_argument("-gm", dest="gmapFiles",
                   nargs="+",
                   required=True,
                   help="Input one or more GMAP (GFF3) files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for TSV-formatted results")
    # Optional - BINge
    p.add_argument("-ga", "--annot", dest="annotationFiles",
                   nargs="+",
                   required=False,
                   help="Optionally, input one or more genome annotation (GFF3) files",
                   default=[])
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many threads to run when multithreading
                   is available (default==1)""",
                   default=1)
    p.add_argument("--convergence_iters", dest="convergenceIters",
                   required=False,
                   type=int,
                   help="""Optionally, specify a maximum number of iterations allowed for
                   bin convergence to be achieved (default==5); in most cases results will
                   converge in fewer than 5 iterations, so setting a maximum acts merely as
                   a safeguard against edge cases I have no reason to believe will ever
                   happen."""
                   if showHiddenArgs else argparse.SUPPRESS,
                   default=5)
    # Optional - program behavioural controls
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
                   (default==0.98)"""
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
    # Help controller
    p.add_argument("--help-long", dest="help-long",
                   action="help",
                   help="""Show all options, including those that are not
                   recommended to be changed."""
                   if not showHiddenArgs else argparse.SUPPRESS)
    
    args = p.parse_args()
    validate_args(args)
    
    # Load indexed transcripts for quick access
    "Load this upfront since it may be a memory limitation, and causing errors early is better than late"
    transcriptRecords = FastaCollection(args.fastaFiles)
    
    # Figure out what the hash of these parameters are, and what our pickle file is called
    paramHash = get_parameters_hash(args)
    pickleFile = os.path.join(
        os.path.dirname(os.path.abspath(args.outputFileName)), 
        f".{paramHash}.binge.pkl"
    )
    
    # Either load a pickle generated by previous BINge run ...
    if os.path.isfile(pickleFile) or os.path.islink(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            binCollection = pickle.load(pickleIn)
    # ... error out if it's a directory or something weird ...
    elif os.path.exists(pickleFile):
        print(f"{pickleFile} already exists, but is not a file?")
        print("BINge expects this to be a file which it can read, or write to.")
        print("Something weird is happening, so I will exit the program now.")
        print(f"Move whatever is at the location of '{pickleFile}' then try again.")
        quit()
    # ... or begin pre-external clustering BINge
    else:
        # Parse each GFF3 into a bin collection structure
        collectionList = generate_bin_collections(args.annotationFiles) # keep genome bins separate to multi-thread later
        
        # Parse GMAP alignments into our bin collection with multiple threads
        novelBinCollection, multiOverlaps = populate_bin_collections(collectionList, args.gmapFiles, args.threads)
        
        # Merge bins resulting from fragmented annotation models
        for i in range(len(collectionList)):
            binCollection = collectionList[i].fix_fragments(multiOverlaps[i])
            collectionList[i] = binCollection
        
        # Merge gene bins together
        binCollection = collectionList[0]
        for i in range(1, len(collectionList)):
            binCollection.merge(collectionList[i])
        
        # Check that this transcriptome file contains the reference gene models
        missingSeqID = find_missing_sequence_id([binCollection, novelBinCollection], transcriptRecords)
        if missingSeqID != None:
            print(f"ERROR: '{missingSeqID}' was binned from GMAP or the annotation GFF3, " +
                "but does not exist in your input transcriptome FASTA.")
            print("A possible error is that your transcriptome lacks the reference sequences.")
            print("This error cannot be reconciled, so the program will exit now.")
            quit()
        
        # Split bins containing overlapping (but not exon-sharing) genes e.g., nested genes
        binCollection = multithread_bin_splitter(binCollection, args.threads)
        
        # Link and merge bins across genomes / across gene copies
        """Usually linking will unify multiple genomes together, but it may detect
        bins of identical gene copies and link them together which is reasonable
        since these would confound DGE to keep separate anyway"""
        
        binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
        novelBinCollection = iterative_bin_self_linking(novelBinCollection, args.convergenceIters)
        
        # Merge bin collections together and re-link bins
        """A novel bin in one genome may just be because of an absence in the annotation.
        Such an absence may not exist in another genome, and hence it will have a gene bin.
        These should be merged together to prevent having redundant bins."""
        
        binCollection.merge(novelBinCollection)
        binCollection = iterative_bin_self_linking(binCollection, args.convergenceIters)
        
        # Write pickle file for potential resuming of program
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(binCollection, pickleOut)
    
    # Write binned clusters to file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#BINge clustering information file\n")
        fileOut.write("cluster_num\tsequence_id\tcluster_type\n")
        
        # Write content lines
        numClusters = 0
        for interval in binCollection:
            bin = interval.data
            for seqID in bin.ids:
                fileOut.write(f"{numClusters+1}\t{seqID}\tbinned\n")
            numClusters += 1
    
    # Cluster remaining unbinned sequences
    unbinnedIDs = get_unbinned_sequence_ids([binCollection], transcriptRecords)
    if len(unbinnedIDs) == 1:
        unbinnedClusterDict = { 0: list(unbinnedIDs)[0] }
    elif len(unbinnedIDs) > 0:
        unbinnedClusterDict = cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args)
    else:
        unbinnedClusterDict = {} # blank to append nothing to output file
    
    # Write output of clustering to file
    with open(args.outputFileName, "a") as fileOut:
        for clusterNum, clusterIDs in unbinnedClusterDict.items():
            for seqID in clusterIDs:
                fileOut.write(f"{clusterNum+numClusters+1}\t{seqID}\tunbinned\n") # clusterType = "unbinned"
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
