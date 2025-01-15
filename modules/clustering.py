import os, sys, random, time
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts.Function_packages import ZS_ClustIO, ZS_BlastIO

def get_hash_for_input_sequences(inObject, randomHash=True, maxLength=20):
    '''
    The function can be guaranteed to produce a truly random hash, or produce something
    consistent. Randomness is good for truly temporary files; consistency can be useful
    for dealing with files that might persist.
    
    Function has been copied from the original in ZS_SeqIO.Conversion to avoid importing
    ZS_SeqIO just for this function.
    
    Parameters:
        inObject -- a FASTA, FastASeq, SeqRecord, or string pointing to an existing FASTA format file
        randomHash -- a boolean indicating whether we want the hash to be consistent when
                    qt is the same (True), or produce truly randomised results always (False)
        maxLength -- an integer for the maximum length of the hash string you want returned
    '''
    assert isinstance(randomHash, bool), \
        "randomHash value must be True or False"
    assert isinstance(maxLength, int) and maxLength > 0, \
        "maxLength must be an integer greater than zero"
    
    # Get a string for hash building
    if isinstance(inObject, str):
        strForHash = inObject
    elif type(inObject).__name__ == "SeqRecord":
        strForHash = inObject.id + str(inObject.seq)[0:1000]
    elif hasattr(inObject, "isFASTA") and inObject.isFASTA is True:
        strForHash = ""
        for FastASeq_obj in inObject:
            strForHash += FastASeq_obj.id
            strForHash += FastASeq_obj.seq[0:1000]
            strForHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    elif hasattr(inObject, "isFastASeq") and inObject.isFastASeq is True:
        strForHash = inObject.id + inObject.seq[0:1000]
    else:
        raise ValueError("get_hash_for_input_sequences can't handle the given object type")
    
    # Get the hash in a randomised or non-randomised way
    if randomHash is True:
        tmpHash = sha256(bytes(strForHash + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
    else:
        tmpHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
    
    return tmpHash[0:maxLength]

def write_unbinned_fasta(unbinnedIDs, transcriptRecords, tmpDir):
    '''
    A helper function which creates a temporary FASTA file containing all unbinned
    sequences. This file will be used for clustering, after which it can be deleted.
    
    Parameters:
        unbinnedIDs -- a set containg sequence IDs that exist in transcriptRecords
                       which should be clustered with an external algorithm.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences.
        tmpDir -- a string location for where this script should write temp files.
    '''
    # Generate a temporary FASTA file containing unbinned transcripts
    tmpFileName = os.path.join(tmpDir, "tmp_BINge_unbinned_{0}.fasta".format(
        get_hash_for_input_sequences(str(transcriptRecords))
        )
    )
    with open(tmpFileName, "w") as fileOut:
        for seqID in unbinnedIDs:
            record = transcriptRecords[seqID]
            fileOut.write(f">{record.name}\n{str(record)}\n")
    
    return tmpFileName

def cluster_unbinned_sequences(unbinnedIDs, transcriptRecords, args, tmpDir, bingeTmpDir):
    '''
    Runs MMSeqs2 or CD-HIT on the unbinned sequences in order to assign them to a cluster.
    It's not ideal, but the alternative is to exclude these sequences which may
    not be in line with the user's aims. By noting the source of clustering in
    the output file, they can decide if they'd like to exclude these or not.
    
    Parameters:
        unbinnedIDs -- a set containg sequence IDs that exist in transcriptRecords
                       which should be clustered with an external algorithm.
        transcriptRecords -- a FASTA file loaded in with pyfaidx for instant lookup of
                             sequences.
        args -- an argparse ArgumentParser object with attributes as set by BINge's
                main argument parsing process.
        tmpDir -- a string location for where MMseqs2 should write temp files.
        bingeTmpDir -- a string location for where the main BINge run is happening
                       and its temporary files should be written.
    '''
    # Generate a temporary FASTA file containing unbinned transcripts
    tmpFileName = write_unbinned_fasta(unbinnedIDs, transcriptRecords, bingeTmpDir)
    
    # Cluster the unbinned transcripts depending on BINge parameters
    if args.unbinnedClusterer in ["mmseqs-cascade", "mmseqs-linclust"]:
        resultClusters = mmseqs_clustering(tmpFileName, args.unbinnedClusterer, 
                                           args.mmseqsDir, tmpDir,
                                           args.threads, args.mmseqsEvalue, args.identity,
                                           args.mmseqsCoverage, args.mmseqsMode,
                                           args.mmseqsSensitivity, args.mmseqsSteps)
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
    mmDB = ZS_BlastIO.MM_DB(fastaFile, mmseqsDir, tmpDir, "nucleotide", threads)
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
    tmpFileName = "tmp_BINge_mms2clusttable_{0}.tsv".format(get_hash_for_input_sequences(fastaFile))
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
