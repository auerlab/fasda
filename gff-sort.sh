#!/bin/sh -e

awk '$3 == "gene" && $1 != "MT" && $1 != "X" && $1 != "Y"' Data/test.gff3 \
    | sort -n -k 1 -k 4 > Data/genes-sorted.gff3
awk '$3 == "mRNA" && $1 != "MT" && $1 != "X" && $1 != "Y"' Data/test.gff3 \
    | sort -n -k 1 -k 4 > Data/mRNA-sorted.gff3
awk '{ print $1 }' Data/test-sorted.gff3 | uniq

