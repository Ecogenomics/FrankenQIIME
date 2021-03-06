#!/usr/bin/env python
# File created on 30 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from optparse import make_option
from glob import glob
from os.path import split
from cogent.util.misc import create_dir
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.format import illumina_data_to_fastq

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Given a directory of per-swath qseq files, this script generates a single fastq per lane."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Generate fastq files from all lanes of read 1 data in the current directory.","process_qseq.py -i ./ -o ./fastq/ -r 1"),
                               ("","Generate fastq files from all lanes of read 2 data in the current directory, truncating the sequences after the first 12 bases.","process_qseq.py -i ./ -o ./fastq/ -r 2 -b 12")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_dir',help='the input directory'),
 make_option('-o','--output_dir',help='the output directory'),
 make_option('-r','--read',help='the read number to consider',type='int')
]
script_info['optional_options'] = [
 make_option('-l','--lanes',
  help='the lane numbers to consider, comma-separated [defaut: %default]',
  default='1,2,3,4,5,6,7,8'),
 make_option('-b','--bases',type='int',
  help='the number of bases to include (useful for slicing a barcode) [defaut: all]',
  default=None),
]
script_info['version'] = __version__

def iter_split_lines(lines):
    """ """
    for line in lines:
        yield tuple(line.split('\t'))

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    input_dir = opts.input_dir
    output_dir = opts.output_dir
    create_dir(output_dir)
    lanes = opts.lanes.split(',')
    bases = opts.bases
    read = opts.read
    
    for lane in lanes:
        read1_fps =  glob('%s/s_%s_%d_*qseq.txt' % (input_dir,
                                                   lane.replace(',',''),
                                                   read))
        # sort so results will be consistent across different runs (important
        # so amplicon and barcodes read headers will match)
        read1_fps.sort()
        for read1_fp in read1_fps:                
            output_fp =  '%s/s_%s_%s_sequences.fastq' % (output_dir,lane,read)
            output_f = open(output_fp,'w')
            for record in iter_split_lines(open(read1_fp,'U')):
                fastq_s = illumina_data_to_fastq(record,
                                                 number_of_bases=bases)
                output_f.write('%s\n' % fastq_s)
            output_f.close()

if __name__ == "__main__":
    main()