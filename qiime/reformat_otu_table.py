#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Adam Skarshewski"
__copyright__ = "Copyright 2012, FrankenQIIME;" 
__credits__ = ["Adam Skarshewski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Adam Skarshewski"
__status__ = "Release"

"""Reformats an OTU table from the QIIME format to the PyroTagger format.
Requires that the FrankenQIIME version of assign_taxonomy.py has been run.
"""

def reformat_otu_table(otu_table_fp, taxonomy_fp):
    
    otu_to_tax = taxonomy_fp_to_dict(taxonomy_fp)
    return write_reformatted_otu_table(otu_table_fp, otu_to_tax)
        
def taxonomy_fp_to_dict(taxonomy_fp):
    
    tax_file = open(taxonomy_fp, 'U')
    otu_to_tax = {}
    
    while True:
        line = tax_file.readline()
        if not line:
            break
        line = line.rstrip()
        if line[0] == '#':
            continue
        splitline = line.split('\t');
        otu_id =  splitline[0]
        taxonomy_stats = {'taxonomy': splitline[1],
                          'evalue': splitline[2],
                          'hit' : splitline[3],
                          'percent_identity': splitline[4],
                          'alignment_length': splitline[5],
                          'mismatches': splitline[6],
                          'gaps': splitline[7],
                          'query':   {'start' : splitline[8],
                                      'end' : splitline[9]},
                          'subject': {'start' : splitline[10],
                                      'end' : splitline[11]}
                           }
        otu_to_tax[otu_id] = taxonomy_stats
    
    tax_file.close()
    return otu_to_tax
    
def write_reformatted_otu_table(otu_table_fp, otu_to_tax):
       
    otu_file = open(otu_table_fp, 'U')
    
    otu_table_string = '#FrankenQIIME ' +  __version__  + " OTU Table\n"
    
    while True:
        line = otu_file.readline()
        if not line:
            break
        line = line.rstrip()
        if line[0] == '#': 
            if line[0:7] == "#OTU ID":
                split_headers = line.split('\t')
                otu_table_string += '\t'.join(split_headers[:-1]) + '\t'
                otu_table_string += '\t'.join(["% Identity", 
                                               "Alignment Length", 
                                               "Mismatches", 
                                               "Gaps", 
                                               "Query Start", 
                                               "Query End", 
                                               "Subject Start", 
                                               "Subject End", 
                                               "E-Value", 
                                               "Hit ID",
                                               "Consensus Lineage"]) + '\n'
            else:
                otu_table_string += line + "\n"
            continue
    
        splitline = line.split('\t');
        otu_id =  splitline[0]
        counts = splitline[1:-1]
        taxonomy_str = splitline[-1]
        
        blast_stats = otu_to_tax[otu_id]
        
        output_stats = [blast_stats['percent_identity'],
                        blast_stats['alignment_length'],
                        blast_stats['mismatches'],
                        blast_stats['gaps'],
                        blast_stats['query']['start'],
                        blast_stats['query']['end'],
                        blast_stats['subject']['start'],
                        blast_stats['subject']['end'],
                        blast_stats['evalue'],
                        blast_stats['hit']
                        ]
        
        otu_table_string += '\t'.join([otu_id] + counts + output_stats + [taxonomy_str]) + '\n'
    
    return otu_table_string