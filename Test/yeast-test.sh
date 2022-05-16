#!/bin/sh -e

cd ..
./cave-man-install.sh
cd Test
dir=Data/05-kallisto-quant
for condition in WT SNF2; do
    ../normalize --output $condition-all-norm.tsv $dir/$condition-*/abundance.tsv
done

