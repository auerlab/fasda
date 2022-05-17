#!/usr/bin/env Rscript

##########################################################################
#   Script description:
#       R script to validate wilcox computations in fold-change.c
#       Hard-coded vectors should match the last output of yeast-test.sh
#       
#   History:
#   Date        Name        Modification
#   2022-05-16  Jason Bacon Begin
##########################################################################

counts1 <- c(6.494263, 15.284300, 3.429830, 11.306324, 13.713963, 16.798491,
	     13.672567, 12.304265, 19.918194, 6.068784, 53.912994, 21.302616,
	     19.949249, 15.349930)
counts1

counts2 <- c(35.113632, 9.991657, 29.798967, 15.793796, 22.247780, 26.676022,
	     25.555389, 21.871218, 40.164658, 14.343984, 53.393482, 53.919132,
	     27.935381, 23.586054)
counts2

wilcox.test(counts1, counts2, alternative="two.sided", exact=TRUE)

U = 33
mU = length(counts1)*length(counts2)/2
sdU = sqrt(length(counts1)*length(counts2)*(length(counts1)+length(counts2)+1)/12)
z = (U-mU)/sdU
print("z:")
print(z)
pval = 2 * pnorm(z)
print("2 * pnorm(z):")
print(pval)

