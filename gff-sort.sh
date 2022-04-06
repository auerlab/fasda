#!/bin/sh -e

awk '$3 == "gene" && $1 != "MT"' Data/test.gff3 \
    | gsort -n -k 1 -k 4 > Data/test-sorted.gff3
awk '{ print $1, $4 }' Data/test-sorted.gff3 | more

