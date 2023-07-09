#!/bin/sh

stringtie -e \
    -G Results/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
    Results/09-hisat-align/SNF2-1.bam \
    -A stringtie.out \
    -o stringtie.gtf

# head 
for transcript in $(cat stringtie.out | cut -f 1); do
    echo $transcript
    awk -v t=$transcript '$1 ~ t { print $4 }' Results/06-kallisto-quant/WT-1/abundance.tsv
    awk -v t=$transcript '$1 == t { print $8 }' stringtie.out
    echo '==='
done | more
