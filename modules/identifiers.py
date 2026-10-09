import os

'''
BINge identifies every sequence by a "sequence key", which is a tuple of
(prefix, seqID). The prefix is the BINge-internal name of the input file that
the sequence came from (e.g., 'genome1', 'annotations2', 'transcriptome1'), and
the seqID is the sequence's identifier within that file. Sequence IDs are only
required to be unique within each input file, so the combination of the two
values is what uniquely identifies a sequence.

Where a sequence must be identified by a single string (e.g., FASTA files given
to external programs like salmon or MMseqs2), the key is flattened into a
namespaced string like 'transcriptome1::seqID'. Since prefixes never contain the
separator, this can always be split back into its original key.
'''
NAMESPACE_SEP = "::"

def prefix_from_file(fileName):
    '''
    Derives the BINge prefix of a sequence or annotation file within the working directory.
    All such files are named like '{prefix}.{suffix}' e.g., 'transcriptome1.cds'.

    Parameters:
        fileName -- a string indicating the location of a file within the BINge working
                    directory.
    Returns:
        prefix -- a string of the file's prefix e.g., 'genome1'.
    '''
    return os.path.basename(fileName).split(".")[0]

def key_to_str(seqKey):
    '''
    Parameters:
        seqKey -- a tuple of (prefix, seqID).
    Returns:
        namespacedID -- a string like '{prefix}::{seqID}'.
    '''
    prefix, seqID = seqKey
    return f"{prefix}{NAMESPACE_SEP}{seqID}"

def str_to_key(namespacedID):
    '''
    Parameters:
        namespacedID -- a string like '{prefix}::{seqID}' as produced by key_to_str().
    Returns:
        seqKey -- a tuple of (prefix, seqID).
    '''
    if not NAMESPACE_SEP in namespacedID:
        raise ValueError(f"'{namespacedID}' is not a namespaced BINge sequence ID (i.e., " +
                         f"'prefix{NAMESPACE_SEP}seqID'); this file may have been generated " +
                         "by an older version of BINge and should be regenerated.")
    prefix, seqID = namespacedID.split(NAMESPACE_SEP, maxsplit=1)
    return (prefix, seqID)
