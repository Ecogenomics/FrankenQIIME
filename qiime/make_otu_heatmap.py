#!/usr/bin/env python
#file make_otu_heatmap.py

from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Release"


from numpy import array,concatenate,asarray,transpose,log,invert,asarray,\
    float32,float64, unique
from cogent.parse.table import SeparatorFormatParser
from optparse import OptionParser
from qiime.parse import parse_otu_table
from qiime.util import MissingFileError
import os
from matplotlib import use
use('Agg',warn=False)
import matplotlib
from matplotlib.pylab import *
from qiime.beta_diversity import get_nonphylogenetic_metric
from cogent.core.tree import PhyloNode
from cogent.cluster.UPGMA import UPGMA_cluster
from qiime.parse import parse_newick, PhyloNode


def get_overlapping_samples(otu_sample_ids, map_rows, otu_table):
    """Extracts only samples contained in otu table and mapping file.
    
       Returns: new_sample_ids, new_map_rows, new_otu_table
    """
    map_sample_ids = zip(*map_rows)[0]
    new_otu_cols = []
    new_map = []
    new_sample_ids = []
    for i, _id in enumerate(otu_sample_ids):
        if _id in map_sample_ids:
            ix = map_sample_ids.index(_id)
            new_otu_cols.append(otu_table[:,i])
            new_sample_ids.append(_id)
            new_map.append(map_rows[ix])
    return array(new_sample_ids), new_map, array(new_otu_cols).T

def extract_metadata_column(sample_ids, metadata, category):
    """Extracts values from the given metadata column"""
    col_ix = metadata[1].index(category)
    map_sample_ids = zip(*metadata[0])[0]
    category_labels = []
    
    for i,sample_id in enumerate(sample_ids):
        if sample_id in map_sample_ids:
            row_ix = map_sample_ids.index(sample_id)
            entry = metadata[0][row_ix][col_ix]
            category_labels.append(entry)
    return category_labels

def get_order_from_categories(otus, category_labels):
    """Groups samples by category values; clusters within each group"""
    category_labels = array(category_labels)
    sample_order = []

    for label in unique(category_labels):
        label_ix = category_labels==label
        label_ix_ix = get_clusters(otus[:,label_ix], axis='column')
        sample_order += list(nonzero(label_ix)[0][array(label_ix_ix)])
    return array(sample_order)


def get_order_from_tree(ids, tree_text):
    """Returns the indices that would sort ids by tree tip order"""
    tree = parse_newick(tree_text, PhyloNode) 
    ordered_ids = []
    for tip in tree.iterTips():
        if tip.Name in ids:
            ordered_ids.append(tip.Name)
    return names_to_indices(ids, ordered_ids)

def make_otu_labels(otu_ids, lineages, n_levels=1):
    """Returns 'pretty' OTU labels: 'Lineage substring (OTU ID)'
        
       Lineage substring includes the last n_levels lineage levels 
    """
    if len(lineages[0]) > 0:
        otu_labels = []
        for i, lineage in enumerate(lineages):
            if n_levels > len(lineage):
                otu_label = '%s (%s)' %(';'.join(lineage),otu_ids[i])
            else:
                otu_label = '%s (%s)' \
                    %(';'.join(lineage[-n_levels:]), otu_ids[i])
            otu_labels.append(otu_label)
        otu_labels = [lab.replace('"','') for lab in otu_labels]
    else:
        otu_labels = otu_ids
    return otu_labels

def names_to_indices(names, ordered_names):
    """Returns the indices that would sort 'names' like 'ordered_names'
    """
    indices = []
    names_list = list(names)
    for ordered_name in ordered_names:
        if ordered_name in names_list:
            indices.append(names_list.index(ordered_name))
    return array(indices)
    

def get_log_transform(data, eps=None):
    """Returns log10 of the data, setting zero values to eps.
    
       If eps is None, eps is set to 1/2 the smallest nonzero value.
    """
    # ensure data are floats
    data = asarray(data,dtype=float64)

    # set all zero entries to a small value
    if eps is None:
        eps = (data[data>0]).min()/2
    data[data==0] = eps
    data = log10(data)
    return data
    

def get_clusters(x_original, axis=['row','column'][0]):
    """Performs UPGMA clustering using euclidean distances"""
    x = x_original.copy()
    if axis=='column':
        x = x.T
    nr = x.shape[0]
    metric_f = get_nonphylogenetic_metric('euclidean')
    row_dissims = metric_f(x)
    # do upgma - rows
    BIG = 1e305
    row_nodes = map(PhyloNode, map(str,range(nr)))
    for i in range(len(row_dissims)):
        row_dissims[i,i] = BIG
    row_tree = UPGMA_cluster(row_dissims, row_nodes, BIG)
    row_order = [int(tip.Name) for tip in row_tree.iterTips()]
    return row_order


def get_fontsize(numrows):
    """Returns the fontsize needed to make text fit within each row.
    """
    thresholds = [25, 50, 75, 100, 125]
    sizes =      [ 5,  4,  3,   2, 1.5,   1]
    i = 0
    while numrows > thresholds[i]:
        i += 1
        if i == len(thresholds):
            break
    return sizes[i]


def plot_heatmap(x, row_labels, col_labels, filename='heatmap.pdf',
        width=5, height=5, textborder=.25):
    """Create a heatmap plot, save as a pdf.
    
        'width', 'height' are in inches
        
        'textborder' is the fraction of the figure allocated for the
        tick labels on the x and y axes
    """
    nrow = x.shape[0]
    ncol = x.shape[1]

    # determine appropriate font sizes for tick labels
    row_fontsize = get_fontsize(nrow)
    col_fontsize = get_fontsize(ncol)

    # create figure and plot heatmap
    fig = figure(figsize=(width, height))
    my_cmap=get_cmap('gist_gray')
    imshow(x[:,::-1],interpolation='nearest', aspect='auto', cmap=my_cmap)
    ax = fig.axes[0]

    # imshow is offset by .5 for some reason
    xlim(-.5, ncol-.5)
    ylim(-.5, nrow-.5)
    
    # add ticklabels to axes
    xticks(arange(ncol), col_labels[::-1], fontsize=col_fontsize)
    yticks(arange(nrow), row_labels, fontsize=row_fontsize)

    # turn off tick marks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # rotate x ticklabels
    for label in ax.xaxis.get_ticklabels():
        label.set_rotation(90)

    # add space for tick labels
    fig.subplots_adjust(left=textborder, bottom=textborder)
    cb = colorbar() # grab the Colorbar instance
    # set colorbar tick labels to a reasonable value (normal is large)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(5)
    fig.savefig(filename)
