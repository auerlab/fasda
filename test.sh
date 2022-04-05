#!/bin/sh -e

make
./diffanal \
    Data/test.gff3 \
    Data/chondro-sample1-rep1-time1/pseudoalignments.bam \
    Data/chondro-sample2-rep1-time2/pseudoalignments.bam
