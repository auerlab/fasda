#!/bin/sh -e

cut -f 1 Results/06-kallisto-quant/cond1-rep01/abundance.tsv > counts.tsv
sample=1
for file in Results/06-kallisto-quant/*/abundance.tsv; do
    cut -f 4 $file | sed -e "s|est_counts|s$sample|" | paste counts.tsv - > temp.tsv
    mv temp.tsv counts.tsv
    sample=$(($sample + 1))
done
head counts.tsv
./normalize.R
