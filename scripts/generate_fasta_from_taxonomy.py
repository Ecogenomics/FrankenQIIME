#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Adam Skarshewski"
__copyright__ = "Copyright 2011, The QIIME Project; 2012, FrankenQIIME"
__credits__ = ["Rob Knight", "Justin Kuczynski","Jesse Stombaugh", "Adam Skarshewski"]
__license__ = "GPL"
__version__ = "1.1.0"
__status__ = "Release"
 
from sys import stdout
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.generate_fasta_from_taxonomy import generate_fasta

options_lookup = get_options_lookup()


script_info={}
script_info['brief_description']="""Generate best BLAST FASTA"""
script_info['script_description']="""The script generate_best_blast_fasta.py takes a taxonomic mapping file (produced by assign_taxonomy.py) and a reference FASTA file to produce a fasta file of best blast hits to the OTU in the taxonomic mapping file."""
script_info['script_usage']=[]
script_info['script_usage'].append(("","""Make an OTU table from an OTU map (i.e., result from pick_otus.py) and a taxonomy assignment file (i.e., result from assign_taxonomy.py). Write the output file to otu_table.txt.""","""%prog -i otu_map.txt -t tax_assignments.txt -o otu_table.txt"""))
script_info['script_usage'].append(("","""Make an OTU table, excluding the sequences listed in chimeric_seqs.txt""","%prog -i otu_map.txt -o otu_table.txt -e chimeric_seqs.txt"))
script_info['output_description']="""The output of make_otu_table.py is a tab-delimited text file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU."""
script_info['required_options']=[\
    make_option('-t','--taxonomy_fp',type="existing_filepath",
    help='path to the input taxonomy file'),
    make_option('-r','--rep_set_fasta_fp',type="existing_filepath",
    help='path to the input fasta file'),
    make_option('-b','--blast_fasta_fp',type="existing_filepath",
    help='path to the input fasta file used to BLAST')
]
script_info['optional_options']=[ \
  options_lookup['output_fp']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    if opts.output_fp:
        outfile = open(opts.output_fp, 'w')
    else:
        outfile = stdout

    outfile.write(generate_fasta(opts.taxonomy_fp, opts.rep_set_fasta_fp, opts.blast_fasta_fp))
    

if __name__ == "__main__":
    main()
