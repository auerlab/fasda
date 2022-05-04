#!/bin/sh -e

bams="Data/Hisat2/chondro-sample1-rep1-time1.bam \
    Data/Hisat2/chondro-sample2-rep1-time2.bam \
    Data/Hisat2/chondro-sample3-rep1-time3.bam"

./cave-man-install.sh

for bam in $bams; do
    printf "Processing $bam...\n"
    base=$(basename $bam)
    raw_abundance=${base%.bam}-abundance.tsv
    echo $raw_abundance
    test -e $raw_abundance || time ./abundance "$@" Data/mRNA-sorted.gff3 $bam
    
    norm_tsv=${base%.bam}-norm.tsv
    time ./normalize < $raw_abundance > $norm_tsv
    
    norm_tsvs="$norm_tsvs $norm_tsv"
done

# chondro-sample1-rep1-time1-norm.tsv
time ./fold-change $norm_tsvs
