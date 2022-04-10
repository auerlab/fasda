#!/bin/sh -e

make
time ./diffanal \
    Data/test-sorted.gff3 \
    Data/Hisat2/chondro-sample1-rep1-time1.bam \
    Data/Hisat2/chondro-sample2-rep1-time2.bam \
    Data/Hisat2/chondro-sample3-rep1-time3.bam \
    2>&1 | tee fold-change.txt
