#!/bin/sh -e

cut -f 1 Results/10-kallisto-quant/sample01*/abundance.tsv > counts.tsv
sample=1
for file in Results/10-kallisto-quant/sample*/abundance.tsv; do
    cut -f 4 $file | sed -e "s|est_counts|s$sample|" | paste counts.tsv - > temp.tsv
    mv temp.tsv counts.tsv
    sample=$(($sample + 1))
done
head counts.tsv
./normalize.R|more
