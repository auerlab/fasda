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

print("counts1:")
counts1 <- as.numeric(read.delim("c1.tsv", header=FALSE))
counts1

print("counts2:")
counts2 <- as.numeric(read.delim("c2.tsv", header=FALSE))
counts2

wilcox.test(counts1, counts2, alternative="two.sided",
	    exact=FALSE, correct=FALSE)

quit()

U = 33
mU = length(counts1)*length(counts2)/2
sdU = sqrt(length(counts1)*length(counts2)*(length(counts1)+length(counts2)+1)/12)
z = (U-mU)/sdU
print("z:")
print(z)
pval = 2 * pnorm(z)
print("2 * pnorm(z):")
print(pval)

