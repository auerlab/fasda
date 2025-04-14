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

# row.names=1 removes "target_id" label from col 1, so it is not
# counted as a data column.
raw_counts = read.delim("counts.tsv", row.names=1)
print("Raw counts:")
head(raw_counts)

# Note: Does not include feature names shown as 1st col
cols = ncol(raw_counts)
print(c("cols = ", cols))

print("Log counts:")
log_counts = log(raw_counts)
head(log_counts)

# Find the average log(count) for each feature (row) and add it to the
# matrix as another column called "mean".  Then we can remove
# entire rows (features) with -Inf in the pseudo_ref column.
#             tibble::rownames_to_column('target_id') %>% 
log_counts = log_counts %>% mutate (mean = rowMeans(log_counts))
print("Log counts with means:")
head(log_counts)
nrow(log_counts)

# Remove features (rows) with -Inf average (pseudo_reference)
filtered_log_counts = log_counts %>% filter(mean != "-Inf")
print("Filtered log counts:")
head(filtered_log_counts)
nrow(filtered_log_counts)

# Subtract avg for feature from log counts
# Note that subtracting from log(x) is the same as dividing x, hence
# we call these ratios, not differences
# 1 means operate on rows
ratio_counts = sweep(filtered_log_counts[,1:cols], 1,
		     filtered_log_counts[,cols+1], "-")
print("Ratio counts (count - mean for this row):")
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
manually_normalized = sweep(raw_counts, 2, scaling_factors, "/")
head(manually_normalized)

mrn = function(raw_counts)
{
    print("Running mrn()...")
    print("Raw counts")
    print(head(raw_counts))
    
    # take the log
    log_counts = log(raw_counts)
    print("Log counts:")
    print(head(log_counts))
    
    # find the psuedo-references per sample by taking the geometric mean
    #             rownames_to_column('target_id') %>% 
    log_counts = log_counts %>% 
		 mutate (mean = rowMeans(log_counts))
    print("Log counts with means:")
    print(head(log_counts))
    
    # Remove features (rows) with -Inf average (pseudo_reference)
    log_counts = log_counts %>% filter(mean != "-Inf")
    print(head(log_counts))
    print(ncol(log_counts))
    print(nrow(log_counts))
    
    # the last columns is the pseudo-reference column 
    mean_column = ncol(log_counts)
    print(c('mean_column = ', mean_column))
    
    # where to stop before the mean column 
    before_mean = mean_column - 1
    print(c('before_mean = ', before_mean))
    
    print(head(log_counts[,1:before_mean]))
    print(head(log_counts[,mean_column]))
    # find the ratio of the log counts to the pseudo-reference
    ratios = sweep(log_counts[,1:before_mean], 1,
		   log_counts[,mean_column], "-")
    print("Ratios:")
    print(head(ratios))
    
    # find the median of the ratios
    sample_medians = apply(ratios, 2, median)
    print("Sample medians:")
    print(head(sample_medians))
    
    # convert the median to a scaling factor
    scaling_factors = exp(sample_medians)
    print("Scaling factors:")
    print(head(scaling_factors))
    
    # use scaling factors to scale the original data
    normalized = sweep(raw_counts, 2, scaling_factors, "/")
    return(normalized)
}

# Compare function results to main
function_normalized = mrn(raw_counts)
head(function_normalized == manually_normalized)
tail(function_normalized == manually_normalized)

print("Normalizing with DESeq2...")
library(DESeq2)

# samples (columns names) of the data should be named
samples = as.data.frame(colnames(raw_counts))
print(samples)

# create a DESeqDataSet object. The design can be altered based on experimental design. A design of 1 means no design. 
head(raw_counts)
dds = DESeqDataSetFromMatrix(countData = raw_counts, colData = samples, design = ~1)
print(head(dds))
exit

# this function generates the size factors
dds = estimateSizeFactors(dds)

# scaling_factors were manually computed using our mor_normalization function
# sizeFactors(dds) is used to find the scaling factors from DESeq2
scaling_factors == sizeFactors(dds)

normalized_deseq2 = counts(dds, normalized = TRUE)

head(normalized_deseq2 == manually_normalized)
tail(normalized_deseq2 == manually_normalized)
