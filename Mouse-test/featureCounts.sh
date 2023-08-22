#!/bin/sh -e

featureCounts -p \
    -a Results/04-reference/Mus_musculus.GRCm39.110.chr.gtf \
    -o counts.txt \
    Results/09-hisat2-align/cond1-rep1.bam

