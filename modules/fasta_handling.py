import os, sys
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts import ZS_SeqIO # expose this to any callers

class FastaCollection:
    '''
    Wrapper for pyfaidx Fasta objects which allows multiple to be combined
    and queried as one logical entity.
    
    Parameters:
        fastaFiles -- a list of strings pointing to the locations of FASTA files
                      which are to be loaded in using pyfaidx.Fasta
    '''
    def __init__(self, fastaFiles):
        self.fastaFiles = fastaFiles
        self.records = []
        
        self._parse_fastas()
    
    def _parse_fastas(self):
        for fastaFile in self.fastaFiles:
            self.records.append(Fasta(fastaFile))
    
    def __getitem__(self, key):
        for records in self.records:
            try:
                return records[key]
            except:
                pass
        raise KeyError(f"'{key}' not found in collection")

    def __iter__(self):
        for records in self.records:
            yield from records
    
    def __repr__(self):
        return (f"<FastaCollection object;num_records='{len(self.records)}';" +
                f"fastaFiles={self.fastaFiles}"
        )
