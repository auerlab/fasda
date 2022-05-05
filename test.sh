#!/bin/sh -e

data_dir=Data/Hisat2
bams="$data_dir/chondro-sample1-rep1-time1.bam \
    $data_dir/chondro-sample2-rep1-time2.bam \
    $data_dir/chondro-sample3-rep1-time3.bam"

./cave-man-install.sh

for bam in $bams; do
    printf "Processing $bam...\n"
    raw_abundance=${bam%.bam}-abundance.tsv
    echo $raw_abundance
    test -e $raw_abundance || time ./abundance "$@" Data/mRNA-sorted.gff3 $bam
    
    norm_tsv=${bam%.bam}-norm.tsv
    time ./normalize < $raw_abundance > $norm_tsv
    
    norm_tsvs="$norm_tsvs $norm_tsv"
done

# chondro-sample1-rep1-time1-norm.tsv
time ./fold-change $norm_tsvs
