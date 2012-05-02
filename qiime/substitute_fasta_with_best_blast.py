#!/usr/bin/env python
__author__ = "Adam Skarshewski"
__copyright__ = "Copyright 2011, The QIIME Project; Copyright 2012, FrankenQIIME"
__credits__ = ["Adam Skarshewski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__status__ = "Release"

from qiime.reformat_otu_table import taxonomy_fp_to_dict
from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs, DNA
import subprocess
from tempfile import TemporaryFile, NamedTemporaryFile
from cogent.app.blast import BlastResult
from time import sleep


max_blast_chunk = 1000

def generate_blast_fasta(query_fasta_fp, ref_fasta_fp, params):
    
    if query_fasta_fp:
            # Get a seq iterator
        seqs = MinimalFastaParser(open(query_fasta_fp))
        # Build object to keep track of the current set of sequence to be
        # blasted, and the results (i.e., seq_id -> (taxonomy,quaility score) 
        # mapping)
    current_seqs = []
    result = {}
    
    # Iterate over the (seq_id, seq) pairs
    for seq_id, seq in seqs:
        # append the current seq_id,seq to list of seqs to be blasted
        current_seqs.append((seq_id,seq))
        
        # When there are 1000 in the list, blast them
        if len(current_seqs) == max_blast_chunk:
            # update the result object
            result.update(get_blast(current_seqs, params))
            # reset the list of seqs to be blasted
            current_seqs = []
    # Assign taxonomy to the remaining sequences
    result.update(get_blast(current_seqs, params))
    
    
    
def get_blast(seqs, params):

        blast_hits = get_blast_hits(seqs, params)
        
        print blast_hits

        # select the best blast hit for each query sequence
        best_blast_hit_ids = self._get_first_blast_hit_per_seq(blast_hits)
        
def get_blast_hits(seqs, params):
        """ blast each seq in seqs against blast_db and retain good hits
        """
        max_evalue = 1e-30
        min_percent_identity = 0.90
        seq_ids = [s[0] for s in seqs]
        result = {}
        blast_params = {'-d': params['blast_db'],'-n':'T', '-a': params['Threads']}
        blast_result = megablast_seqs(seqs, blast_params)    
        for seq_id in seq_ids:
            blast_result_id = seq_id.split()[0]
            try:
                result[seq_id] = [
                e for e in blast_result[blast_result_id][0]
                 if (float(e['E-VALUE']) <= max_evalue and \
                  float(e['% IDENTITY']) >= min_percent_identity)]
            except KeyError:
                result[seq_id] = []

        return result
        
def get_first_blast_hit_per_seq(self, blast_hits):
    """ discard all blast hits except the best for each query sequence
    """
    result = {}
    for k,v in blast_hits.items():
        k = k.split()[0]    #get rid of spaces
        try:
            result[k] = v[0]
        except IndexError:
            # If there is no good blast hit, do we want to 
            # leave the key out, or have it point to None?
            result[k] = None
    
    return result



def megablast_seqs(seqs, params):
    seq_temp = NamedTemporaryFile(bufsize=0)
    seq_temp.write(LoadSeqs(data=seqs, moltype=DNA, aligned=False).toFasta())
    seq_temp.flush()
    out_temp = TemporaryFile()
    
    print seq_temp.name    
    
    params['-i'] = seq_temp.name
    params['-m'] = 9
    print params
    param_list = reduce((lambda x, y: x + y),
        [[str(k), str(v)] for (k, v) in params.items()])
    
    print param_list
    proc_handle = subprocess.Popen(["megablast"] + param_list, stdout=out_temp)
    proc_handle.wait()
    
    lines = [line for line in out_temp]
    print lines
    blast_result = BlastResult(lines)
    
    seq_temp.close()
    out_temp.close()
    
    return blast_result
