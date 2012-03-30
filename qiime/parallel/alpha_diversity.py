#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# alpha_diversity.py

from __future__ import division
from os.path import split, splitext
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

def get_job_commands(python_exe_fp,alpha_diversity_fp,tree_fp,job_prefix,\
    metrics,input_fps,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate alpha diversity commands to be submitted to cluster
    """

    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    commands = []
    result_filepaths = []
    
    for input_fp in input_fps:
        input_path, input_fn = split(input_fp)
        output_fn = 'alpha_%s' % input_fn
        rename_command, current_result_filepaths = get_rename_command(\
         [output_fn],working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -i %s -o %s -t %s -m %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          alpha_diversity_fp,\
          input_fp,
          working_dir + '/' + output_fn,
          tree_fp,
          metrics,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths
