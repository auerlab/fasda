#!/usr/bin/env Rscript

##########################################################################
#   Script description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2025-04-13  Jason Bacon Begin
##########################################################################

# library()
# load()
library(dplyr)
library(tibble)

counts = read.delim("counts.tsv", row.names=1)
head(counts)
log_counts = log(counts)
head(log_counts)

# Find the average log(count) for each feature
log_counts = log_counts %>% 
	     rownames_to_column('gene') %>% 
	     mutate (pseudo_reference = rowMeans(log_counts))
head(log_counts)

# Remove features with -Inf average
filtered_log_counts = log_counts %>% filter(pseudo_reference != "-Inf")
head(filtered_log_counts)

# Subtract avg for feature from log counts
# FIXME: Verify that 8 is correct
ratio_counts = sweep(filtered_log_counts[,2:7], 1, filtered_log_counts[,8], "-")
head(ratio_counts)

# Find median of ratios for each sample
sample_medians = apply(ratio_counts, 2, median)
sample_medians

# Transform back from log space to counts to get scaling factor
scaling_factors = exp(sample_medians)
scaling_factors

# Divide original counts by scaling factor for sample
manually_normalized = sweep(counts, 2, scaling_factors, "/")
head(manually_normalized)

