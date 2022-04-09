#!/bin/sh -e

make
./diffanal \
    Data/test-sorted.gff3 \
    Data/chondro-sample1-rep1-time1/pseudoalignments.bam \
    Data/chondro-sample2-rep1-time2/pseudoalignments.bam \
    Data/chondro-sample3-rep1-time3/pseudoalignments.bam
