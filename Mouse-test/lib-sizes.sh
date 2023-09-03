#!/bin/sh -e

for fastq in Results/01-fetch/Raw-renamed/cond*.zst; do
    printf "$fastq: "
    zstdcat $fastq | wc -l
done
