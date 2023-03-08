#!/bin/sh -e

cd Data/06-kallisto-quant
genes=$(head -n 11 SNF2-1/abundance.tsv | tail -10 | cut -f 1)
for gene in $genes; do
    awk -v gene=$gene '$1 == gene' WT-*/abundance.tsv | cut -f 5 | sort -g > temp
    min=$(head -1 temp)
    median=$(head -24 temp | tail -1)
    max=$(tail -1 temp)
    printf "%s\t%s\t%s\t%s\n" $gene $min $median $max
done
