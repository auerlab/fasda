#!/usr/bin/env Rscript

##########################################################################
#   Script description:
#       Perform median of ratios normalization and DESeq2 analysis based on:
#       https://scienceparkstudygroup.github.io/\
#           research-data-management-lesson/\
#           median_of_ratios_manual_normalization/index.html
#       
#   History:
#   Date        Name        Modification
#   2025-04-13  Jason Bacon Begin
##########################################################################

library(dplyr)
library(tibble)

##########################################################################
# Input file must be TSV and look like this:
# target_id       s1      s2      s3      s4      s5      s6
# YPL071C_mRNA    14      32      30      51      36      36
# YLL050C_mRNA    238     433     344     831     676     631
# YMR172W_mRNA    42.6436 62.4375 46.8129 121.99  100.89  103.369
#
# Column labels are not important.
# First 3 columns are condition 1 (control in DESeq2 design)
# Last 3 are condition 2 (treated in DESeq2 design)
##########################################################################

# row.names=1 removes "target_id" label from col 1, so it is not
# counted as a data column.
raw_counts = read.delim("counts.tsv", row.names=1, header=TRUE)
print("Raw counts:")
head(raw_counts)
cols = ncol(raw_counts)
print(paste("cols = ", cols))

# Transform into log() space to simplify calculations and filtering
print("Log counts:")
log_counts = log(raw_counts)
head(log_counts)

# Find the mean log(count) for each feature (row) and append to the
# matrix as another column called "mean".  Then we can remove
# entire rows (features) with -Inf in the mean column.
# Mean is also called the pseudo_reference in some docs.
# Removed: tibble::rownames_to_column('target_id') %>% 
log_counts = log_counts %>% mutate (mean = rowMeans(log_counts))
print("Log counts with means:")
head(log_counts)
nrow(log_counts)

# Remove features (rows) with -Inf mean log(count) (pseudo_reference)
# Same as geometric mean of counts == 0?
filtered_log_counts = log_counts %>% filter(mean != "-Inf")
print("Filtered log counts:")
head(filtered_log_counts)
nrow(filtered_log_counts)

# Subtract mean for feature from all log counts for feature (row).
# Note that subtracting from log(x) is the same as dividing x, hence
# we call these ratios, not differences.
# 1 means operate on rows.
ratio_counts = sweep(filtered_log_counts[,1:cols], 1,
		     filtered_log_counts[,cols+1], "-")
print("Ratios (log(count) - mean log(count) for this row):")
head(ratio_counts)

# Find median of all "ratios" for each sample
# Store in a 1-row data frame
# 2 means operate on columns
# median is an R function
sample_medians = apply(ratio_counts, 2, median)
sample_medians

# Transform back from log space to get raw counts scaling factor
scaling_factors = exp(sample_medians)
scaling_factors

# Divide raw counts by scaling factor for sample
# 2 means operate on columns
normalized = sweep(raw_counts, 2, scaling_factors, "/")
head(normalized)

##########################################################################
#   Function to do the same as above
##########################################################################

mrn = function(raw_counts)
{
    log_counts = log(raw_counts)
    log_counts = log_counts %>% 
		 mutate (mean = rowMeans(log_counts))
    log_counts = log_counts %>% filter(mean != "-Inf")
    mean_column = ncol(log_counts)
    before_mean = mean_column - 1
    ratios = sweep(log_counts[,1:before_mean], 1,
		   log_counts[,mean_column], "-")
    sample_medians = apply(ratios, 2, median)
    scaling_factors = exp(sample_medians)
    normalized = sweep(raw_counts, 2, scaling_factors, "/")
    return(normalized)
}

# Compare function results to results from main
function_normalized = mrn(raw_counts)
head(function_normalized == normalized)
tail(function_normalized == normalized)

##########################################################################
#   Compare DESeq2 normalization
#   https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
##########################################################################

print("Normalizing with DESeq2...")
library(DESeq2)

# DESeqDataSetFromMatrix() can only take integers.  Seriously??
raw_counts = round(raw_counts)

# Create a data frame describing which columns represent each condition.
# For a thorough and comprehensible tutorial:
# https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial
sample_names = c(colnames(raw_counts))
# Must be called "condition" for DESeqDataSetFromMatrix()
condition = c("control", "control", "control", "treated", "treated", "treated")
meta_data = data.frame(sample_names, condition)
# Drop "1", "2", ... and just keep sample names and associated conditions
meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="sample_names")
meta_data

# Sanity checks: Should print TRUE TRUE
all(colnames(raw_counts) %in% rownames(meta_data))
all(colnames(raw_counts) == rownames(meta_data))

# create a DESeqDataSet object.  See docs for possible designs.
# Typical is a 2-condition experiment.
dds = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta_data,
			     design = ~condition)
dds

# Compute normalized counts: Not used by DESeq2, just informational.
# Compute median ratio normalization scaling factors
dds = estimateSizeFactors(dds)
# Extract scaling factors from DESeqDataSet into standard data frame.
scaling_factors == sizeFactors(dds)
scaling_factors
normalized_deseq2 = counts(dds, normalized = TRUE)
write.table(normalized_deseq2, file='deseq2-normalized-counts.tsv',
	    sep = '\t', quote=F, col.names = NA)

# Since we had to round the counts to integers to let DESeq2 use them,
# they won't be quite the same as the manually normalized counts.  Rather
# than compare them, which will yield mostly FALSE values, just display
# them for visual comparison.  They should be about the same.
head(normalized)
head(normalized_deseq2)
tail(normalized)
tail(normalized_deseq2)

# Run DESeq2 differential analysis to produce fold-changes and P-values.
dds = DESeq(dds)
results = results(dds)
results
write.table(results, file='deseq2-results.tsv', sep='\t', quote=F, col.names=NA)
