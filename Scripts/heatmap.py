#!/usr/bin/env python

# Derived from https://igb.mit.edu/mini-courses/python/data-processing-with-python/seaborn/visualizing-rnaseq-data

import sys, os
import pandas
import numpy
import scipy
import seaborn
import glob
import matplotlib.pyplot as plt
import fastcluster

#############################################################################
# Requires one command-line argument containing header line and
# counts following feature name.  Columns should be
# condition1-rep1, condition1-rep2, ... condition2-rep1, ...
#
# Feature c1r1    c1r2    c1r3    c2r1    c2r2    c2r3
# YKL103C_mRNA    157.915324      117.014504      171.523362      558.813838      380.590665      343.314298
# YIR015W_mRNA    11.759652       12.317316       7.146807        21.891676       26.574909       15.902842
#############################################################################

def usage():
    print("Usage: %s %s" % (sys.argv[0], "[--debug] counts-file.tsv"))
    sys.exit(1)

debug = False
if len(sys.argv) == 2:
    counts_file = sys.argv[1]
elif len(sys.argv) == 3:
    if sys.argv[1] == '--debug':
        debug = True
        counts_file = sys.argv[2]
    else:
        usage()
else:
    usage()

seaborn.set_context('paper')
seaborn.set_style("whitegrid")

counts = pandas.read_csv(counts_file,sep='\t')
if debug:
    print(counts.shape)
    print(counts.head)
    print(counts.describe())

row_medians = counts.median(axis=1,numeric_only=True)
if debug:
    print(row_medians)

# Normalize for heatmap
counts_row_norm = counts.copy()
for col in counts.columns[1:]:
    counts_row_norm[col] = numpy.log2((counts[col]+0.1)/(row_medians+0.1))
if debug:
    print(counts_row_norm.describe())
    # print(counts_row_norm.iloc[6,1:].median()) #test some random rows, make sure the median value is 0
    print(counts_row_norm)

counts_row_norm = counts_row_norm.set_index('Feature')
# seaborn.clustermap(counts_row_norm)
seaborn.clustermap(counts_row_norm, col_cluster=False, cmap="RdBu_r", figsize=(10,10))
plt.show()
