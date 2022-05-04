#!/bin/sh -e

./cave-man-install.sh
time ./abundance "$@" \
    Data/mRNA-sorted.gff3 \
    Data/Hisat2/chondro-sample1-rep1-time1.bam \
    Data/Hisat2/chondro-sample2-rep1-time2.bam \
    Data/Hisat2/chondro-sample3-rep1-time3.bam
