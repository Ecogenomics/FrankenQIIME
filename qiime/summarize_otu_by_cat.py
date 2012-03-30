#!/usr/bin/env python
#file summarize_otu_by_cat.py
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Julia Goodrich", "Jesse Stombaugh","Greg Caporaso","Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Release"


"""
Author: Julia Goodrich (julia.goodrich@colorado.edu)
Status: Prototype

Requirements:
Python 2.5

This script generates the otu table for a specific category
"""

from optparse import OptionParser
from collections import defaultdict
from numpy import nonzero, arange, array
from string import strip
import os
from time import strftime
from random import choice, randrange
from qiime.format import format_otu_table
from decimal import getcontext
from qiime.parse import parse_mapping_file, parse_otu_table

def get_sample_cat_info(lines, category):
    cat_by_sample = {}
    sample_by_cat = defaultdict(list)
    meta_dict = {}
    num_samples_by_cat = defaultdict(int)
    label_lists_dict = defaultdict(list)
    mapping_data, header, comments = parse_mapping_file(lines)
    
    category_labels = header
    index = category_labels.index(category)

    for line in mapping_data:
        categories = line[0:len(category_labels)+1]
        sample = categories[0].strip()
        meta_dict[sample] = [(categories[index],0)]

        cat_by_sample[sample] = [(l.strip(),c.strip()) \
                             for l,c in zip(category_labels,categories)]

        cat_list = []
        for i,(l,c) in enumerate(zip(category_labels,categories)):
            if c not in label_lists_dict[l]:
                label_lists_dict[l].append(c)
            l = l.strip()
            c = c.strip()
            cat_list.append((l,c))
            sample_by_cat[(l,c)].append(sample)
            num_samples_by_cat[(l,c)] += 1

        cat_by_sample[sample] = cat_list

    return cat_by_sample, sample_by_cat, len(category_labels), meta_dict,label_lists_dict,num_samples_by_cat


def get_counts_by_cat(lines, num_meta, meta_dict, cat_list,category,num_samples_by_cat,
                     normalize):
    con_by_sample = defaultdict(set)
    node_file_str = []
    edge_file_str = []
    red_nodes = defaultdict(int)
    red_node_file_str = []
    red_edge_file_str = []
    edge_from = []
    to = []
    otu_dc = defaultdict(int)
    degree_counts = defaultdict(int)
    sample_dc = defaultdict(int)
    sample_num_seq = defaultdict(int)
    samples_from_mapping = meta_dict.keys()
    con_list = []
    label_list = []
    norm_otu_table =[]
    sample_counts = defaultdict(int)
    cat_otu_table = []
    otus = []
    taxonomy = []
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(lines)

    label_list = sample_ids
    if lineages == []:
        is_con = False
    else:
        is_con = True
    for idx, line in enumerate(otu_table):
        new_line = []
        label_dict = defaultdict(int)
        data = line
        to_otu = otu_ids[idx]
        otus.append(to_otu)
        con = ''
        if is_con:
            con = '; '.join(lineages[idx])
            counts = data
        else:
            counts = data
        taxonomy.append(con)
        if not normalize:
            for i,c in zip(label_list,counts):
                if i in samples_from_mapping:
                    label_dict[meta_dict[i][0][0]] += c        
            for i in cat_list:
                new_line.append(str(label_dict[i]))
            cat_otu_table.append(new_line)

        else:
            new_line.extend(counts)
            norm_otu_table.append(new_line)
            for i, c in zip(label_list,counts):
                sample_counts[i] += c
    total = 0
    if normalize:
        for l in norm_otu_table:
            counts = l
            new_line = []
            label_dict = defaultdict(float)
            getcontext().prec = 28
            for i,c in zip(label_list,counts):
                if i in samples_from_mapping:
                    label_dict[meta_dict[i][0][0]] += float(c)/(sample_counts[i])
            for i in cat_list:
                new_line.append(round((label_dict[i]/ num_samples_by_cat[(category,i)]),5))
            cat_otu_table.append(new_line)
    return  cat_otu_table, otus, taxonomy


def summarize_by_cat(map_lines,otu_sample_lines,category,norm):
    """creates the category otu table"""
    cat_by_sample, sample_by_cat, num_meta, meta_dict, label_lists_dict, \
                   num_samples_by_cat = get_sample_cat_info(map_lines,category)

    lines, otus, taxonomy = get_counts_by_cat(otu_sample_lines, num_meta, \
                  meta_dict,label_lists_dict[category],category,num_samples_by_cat,\
                  norm)
    
    #This for loop was added to remove columns that sum to 0, since you may 
    #pass a mapping file that has more samples than in the OTU table, hence resulting
    #in columns with no counts
    new_labels=[]
    new_lines=[]
    for i,line in enumerate(zip(*lines)):
        total_col=sum([float(x) for x in line])
        if total_col>0:
            new_lines.append(line)
            new_labels.append(label_lists_dict[category][i])
    new_lines=zip(*new_lines)
    
    lines = format_otu_table(new_labels, otus, array(new_lines), \
                  taxonomy=taxonomy,
                  comment='Category OTU Counts-%s'% category)
    return lines
