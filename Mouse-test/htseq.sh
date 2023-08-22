#!/bin/sh

# Very slow compared to stringtie and featureCounts
htseq-count Results/09-hisat2-align/cond1-rep1.bam \
    Results/04-reference/Mus_musculus.GRCm39.110.chr.gtf \
    > counts-htseq.txt

