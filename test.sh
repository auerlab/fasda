#!/bin/sh -e

bams="Data/Hisat2/chondro-sample1-rep1-time1.bam \
    Data/Hisat2/chondro-sample2-rep1-time2.bam \
    Data/Hisat2/chondro-sample3-rep1-time3.bam"

./cave-man-install.sh

for bam in $bams; do
    printf "Processing $bam...\n"
    base=$(basename $bam)
    intsv=${base%.bam}-abundance.tsv
    echo $intsv
    test -e $intsv || time ./abundance "$@" Data/mRNA-sorted.gff3 $bam
    
    outtsv=${base%.bam}-norm.tsv
    echo $tsv
    time ./normalize < $intsv > $outtsv
done

# time ./fold-change outtsv1 outtsv2 outtsv3
