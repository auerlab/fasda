#!/usr/bin/env Rscript

##########################################################################
#   Script description:
#       Perform median of ratios normalization in R based on:
#       https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
#       
#   History:
#   Date        Name        Modification
#   2025-04-13  Jason Bacon Begin
##########################################################################

library(dplyr)
library(tibble)

counts = read.delim("counts.tsv", row.names=1)
head(counts)
# Note: Does not include feature names shown as 1st col
cols = ncol(counts)

log_counts = log(counts)
head(log_counts)

# Find the average log(count) for each feature (row) and add it to the
# matrix as another column called "pseudo_reference".  Then we can remove
# entire rows (features) with -Inf in the pseudo_ref column.
log_counts = log_counts %>% 
	     tibble::rownames_to_column('target_id') %>% 
	     mutate (pseudo_reference = rowMeans(log_counts))
head(log_counts)
nrow(log_counts)

# Remove features (rows) with -Inf average (pseudo_reference)
filtered_log_counts = log_counts %>% filter(pseudo_reference != "-Inf")
head(filtered_log_counts)
nrow(filtered_log_counts)

# Subtract avg for feature from log counts
# Note that subtracting from log(x) is the same as dividing x, hence
# we call these ratios, not differences
# 1 means operate on rows
ratio_counts = sweep(filtered_log_counts[,1:cols+1], 1,
		     filtered_log_counts[,cols+2], "-")
print("Ratio counts:")
head(ratio_counts)

# Find median of ratios for each sample
# Store in a 1-row data frame
# 2 means operate on columns
# median is an R function
sample_medians = apply(ratio_counts, 2, median)
sample_medians

# Transform back from log space to counts to get scaling factor
scaling_factors = exp(sample_medians)
scaling_factors

# Divide original counts by scaling factor for sample
# 2 means operate on columns
manually_normalized = sweep(counts, 2, scaling_factors, "/")
head(manually_normalized)

mrn = function(counts)
{
  print("Running mrn()...")
  head(counts)
  
  # take the log
  log_counts = log(counts) 
  head(log_counts)
  
  # find the psuedo-references per sample by taking the geometric mean
  log_counts = log_counts %>% 
	       rownames_to_column('target_id') %>% 
	       mutate (gene_averages = rowMeans(log_counts)) %>% 
	       filter(gene_averages != "-Inf")
  head(log_counts)
  return(log_counts)
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_counts)
  print('pseudo_column = ')
  print(pseudo_column)
  
  # where to stop before the pseudo column 
  before_pseudo = pseudo_column - 1
  
  # find the ratio of the log counts to the pseudo-reference
  ratios = sweep(log_counts[,2:before_pseudo], 1,
		 log_counts[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}

# Not yet working
# manually_normalized == mrn(counts)

