#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.alpha_diversity import (single_file_alpha, multiple_file_alpha,
list_known_metrics)
import os

#alpha_diversity.py
script_info={}
script_info['brief_description']="""Calculate alpha diversity on each sample in an otu table, using a variety of alpha diversity metrics"""
script_info['script_description']="""This script calculates alpha diversity, or within-sample diversity, using an otu table. The QIIME pipeline allows users to conveniently calculate more than two dozen different diversity metrics. The full list of available metrics is available by passing the option -s to the script alpha_diversity.py. Every metric has different strengths and limitations - technical discussion of each metric is readily available online and in ecology textbooks, but is beyond the scope of this document."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Single File Alpha Diversity Example:""","""To perform alpha diversity (e.g. chao1) on a single OTU table, where the results are output to "alpha_div.txt", you can use the following command:""","""alpha_diversity.py -i otu_table.txt -m chao1 -o alpha_div.txt"""))
script_info['script_usage'].append(("""""","""Note: Since this is a non-phylogenetic metric, the tree does not need to be supplied.""",""""""))
script_info['script_usage'].append(("""""","""In the case that you would like to perform alpha diversity using a phylogenetic metric (e.g. PD_whole_tree), you can use the following command:""","""alpha_diversity.py -i otu_table.txt -m PD_whole_tree -o alpha_div.txt -t repr_set.tre"""))
script_info['script_usage'].append(("""""","""You can use the following idiom to run multiple metrics at once (comma-separated):""","""alpha_diversity.py -i otu_table.txt -m chao1,PD_whole_tree -o alpha_div.txt -t repr_set.tre"""))
script_info['script_usage'].append(("""Multiple File (batch) Alpha Diversity:""","""To perform alpha diversity on multiple OTU tables (e.g.: rarefied otu tables resulting from multiple_rarefactions.py), specify an input directory instead of a single otu talbe, and an output directory (e.g. "alpha_div_chao1_PD/") as shown by the following command:""","""alpha_diversity.py -i rarefaction_tables/ -m chao1,PD_whole_tree -o alpha_div_chao1_PD/ -t repr_set.tre"""))
script_info['output_description']="""The resulting file(s) is a tab-delimited text file, where the columns correspond to alpha diversity metrics and the rows correspond to samples and their calculated diversity measurements. When a folder is given as input (-i), the script processes every otu table file in the given folder, and creates a corresponding file in the output directory.

Example Output:

====== ======= ============= ================
\      simpson PD_whole_tree observed_species
====== ======= ============= ================
PC.354 0.925   2.83739       16.0
PC.355 0.915   3.06609       14.0
PC.356 0.945   3.10489       19.0
PC.481 0.945   3.65695       19.0
PC.593 0.91    3.3776        15.0
PC.607 0.92    4.13397       16.0
PC.634 0.9     3.71369       14.0
PC.635 0.94    4.20239       18.0
PC.636 0.925   3.78882       16.0
====== ======= ============= ================
"""
script_info['required_options']=[]
script_info['optional_options']=[\
 make_option('-i', '--input_path',
     help='Input OTU table filepath or input directory containing OTU' +\
     ' tables for batch processing. [default: %default]',
     type='existing_path'),
 make_option('-o', '--output_path',
     help='Output distance matrix filepath or output directory to store' +\
     ' distance matrices when batch processing. [default: %default]',
     type='new_path'),
 make_option('-m', '--metrics', default='PD_whole_tree,chao1,observed_species',
     help='Alpha-diversity metric(s) to use. A comma-separated list should' +\
     ' be provided when multiple metrics are specified. [default: %default]'), 
 make_option('-s', '--show_metrics', action='store_true', 
     dest="show_metrics",
     help='Show the available alpha-diversity metrics and exit.'),
 make_option('-t', '--tree_path', default=None,
     help='Input newick tree filepath.' +\
     ' [default: %default; REQUIRED for phylogenetic metrics]',
     type='existing_filepath')
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.show_metrics:
        print("Known metrics are: %s\n" \
              % (', '.join(list_known_metrics()),))
        exit(0)
    almost_required_options = ['input_path','output_path','metrics']
    for option in almost_required_options:
        if getattr(opts,option) == None:
            option_parser.error('Required option --%s omitted.' % option)
    
    if os.path.isdir(opts.input_path):
      multiple_file_alpha(opts.input_path, opts.output_path, opts.metrics, 
        opts.tree_path)
    elif os.path.isfile(opts.input_path):
      try:
          f = open(opts.output_path, 'w')
          f.close()
      except IOError:
          print("ioerror, couldn't create output file")
          exit(1)
      single_file_alpha(opts.input_path, opts.metrics, 
          opts.output_path, opts.tree_path)
    else:
      print("io error, input path not valid. does it exist?")
      exit(1)

if __name__ == "__main__":
    main()