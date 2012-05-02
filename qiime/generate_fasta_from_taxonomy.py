#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Adam Skarshewski"
__copyright__ = "Copyright 2012, FrankenQIIME;" 
__credits__ = ["Adam Skarshewski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__status__ = "Release"

from qiime.reformat_otu_table import taxonomy_fp_to_dict
from cogent.parse.fasta import MinimalFastaParser

def generate_fasta(taxonomy_fp, rep_set_fasta_fp, blast_fasta_fp):
    
    otu_to_tax = taxonomy_fp_to_dict(taxonomy_fp)
           
    rep_set_header_headers = get_rep_set_headers(rep_set_fasta_fp)
    
    return get_fasta_from_tax(blast_fasta_fp, otu_to_tax,
                              rep_set_header_headers)
   
def get_fasta_from_tax(ref_fasta_file_fp, otu_to_tax, orig_headers):
   
    fasta_file = open(ref_fasta_file_fp)
       
    fasta_str = ''
    blast_seqs = dict((name,seq) for (name,seq) in MinimalFastaParser(fasta_file))
    for otu_id in sorted(orig_headers, key=(lambda x: orig_headers[x]["position"])):
        if otu_id in otu_to_tax:
            # Anything that doesn't have a BLAST hit is omitted
            if otu_to_tax[otu_id]["hit"] == "-":
                continue
            fasta_str += ">%s\n%s\n" % (orig_headers[otu_id]["header"],
                                        blast_seqs[otu_to_tax[otu_id]["hit"]])
    fasta_file.close()
    return fasta_str     

def get_rep_set_headers(rep_set_fasta_fp):
   
    fasta_file = open(rep_set_fasta_fp)
       
    headers = {}
    count = 0
    for (name, seq) in MinimalFastaParser(fasta_file):
        count += 1
        otu_id = name.split(" ")[0]
        headers[otu_id] = {"position": count,
                          "header": name}
        
    fasta_file.close()
    return headers
        
        