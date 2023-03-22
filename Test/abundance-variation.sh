#!/bin/sh -e

cd Results/06-kallisto-quant

# First ten transcripts
genes=$(head -n 11 SNF2-1/abundance.tsv | tail -10 | cut -f 1)

# All transcripts
# genes=$(cat SNF2-1/abundance.tsv | cut -f 1)

printf "%s\t%6s\t%6s\t%6s\t%6s\n" "Transcript" "Min" "Median" "Max" "Spread"
for gene in $genes; do
    awk -v gene=$gene '$1 == gene' WT-*/abundance.tsv | cut -f 5 | sort -g > temp
    min=$(head -1 temp)
    median=$(head -24 temp | tail -1)
    max=$(tail -1 temp)
    spread=$(echo "scale=1; $max / $min" | bc -l)
    printf "%s\t%6.1f\t%6.1f\t%6.1f\t%6.1f\n" $gene $min $median $max $spread
done
