#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project; Copyright 2012, FrankenQIIME;" 
__credits__ = ["Rob Knight", "Justin Kuczynski", "Adam Skarshewski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Adam Skarshewski"
__email__ = "a.skarshewski@uq.edu.au"
__status__ = "Release"

"""Reformats an OTU table from the QIIME format to the PyroTagger format.
Requires that the FrankenQIIME version of assign_taxonomy.py has been run.
"""

from sys import argv, exit, stderr, stdout
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.reformat_otu_table import reformat_otu_table

options_lookup = get_options_lookup()

#make_otu_table.py
script_info={}
script_info['brief_description']="""Reformat OTU table"""
script_info['script_description']="""This script takes a QIIME formatted OTU table and converts it to the FrankenQIIME (PyroTagger) Format."""
script_info['script_usage']=[]
script_info['script_usage'].append(("","""Reformat an QIIME OTU table (i.e., result from make_otu_table.py) using a taxonomy assignment file (i.e., result from assign_taxonomy.py)..""","""%prog -i otu_table.txt -t tax_assignments.txt -o extended_otu_table.txt"""))
script_info['output_description']="""The output of %prog is a tab-delimited text file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU."""
script_info['required_options']=[\
 options_lookup['otu_table_as_primary_input'],
  make_option('-t', '--taxonomy', dest='taxonomy_fname', \
              help='Path to taxonomy assignment, containing the assignments of \ taxons to sequences (i.e., resulting txt file from assign_taxonomy.py) \
 [default: %default]', default=None)
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

    outfile.write(reformat_otu_table(opts.otu_table_fp, opts.taxonomy_fname))
    
if __name__ == "__main__":
    main()
