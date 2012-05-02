#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Adam Skarshewski"
__copyright__ = "Copyright 2012, FrankenQIIME"
__credits__ = ["Adam Skarshewski"]
__license__ = "GPL"
__version__ = "1.1.0"
__status__ = "Release"
 
from sys import stdout
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.substitute_fasta_with_best_blast import generate_blast_fasta

options_lookup = get_options_lookup()


script_info={}
script_info['brief_description']="""Substitute FASTA with Best BLAST"""
script_info['script_description'] = """The script substitute_fasta_with_blast.py takes a query fasta file and a reference fasta file, and substitutes the entries in the query file with their best (mega)BLAST hit in the reference."""
script_info['script_usage']=[]
script_info['script_usage'].append(("","""Make an OTU table from an OTU map (i.e., result from pick_otus.py) and a taxonomy assignment file (i.e., result from assign_taxonomy.py). Write the output file to otu_table.txt.""","""%prog -i otu_map.txt -t tax_assignments.txt -o otu_table.txt"""))
script_info['script_usage'].append(("","""Make an OTU table, excluding the sequences listed in chimeric_seqs.txt""","%prog -i otu_map.txt -o otu_table.txt -e chimeric_seqs.txt"))
script_info['output_description']="""The output of substitute_fasta_with_best_blast.py is a FASTA file, where the entries in the query file have been replaced by the corresponding best BLAST hit in the reference FASTA (or omitted if below the e-value threshold)."""
script_info['required_options']=[\
    make_option('-q','--query_fasta_fp',type="existing_filepath",
    help='path to the query FASTA file'),
    make_option('-r','--ref_fasta_fp',type="existing_filepath",
    help='path to the reference FASTA file'),
    make_option('-d','--ref_blastdb_fp',
        help='path to the blastdb created from the reference FASTA file')
]
script_info['optional_options']=[ \
  options_lookup['output_fp'],
  make_option('-a', '--threads',  type='int',
        help='Number of threads to use for blast [default: %default]',
        default=1)
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    params = {}    
    if opts.output_fp:
        outfile = open(opts.output_fp, 'w')
    else:
        outfile = stdout
    
    params['reference_seqs_filepath'] = opts.ref_fasta_fp
    params['blast_db'] = opts.ref_blastdb_fp
    
    if opts.threads:    
        params['Threads'] = opts.threads        

    outfile.write(generate_blast_fasta(opts.query_fasta_fp,
                                       opts.ref_fasta_fp,
                                       params))
    

if __name__ == "__main__":
    main()
