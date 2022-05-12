#!/bin/sh -e

mkdir -p Data/Yeast/Raw-renamed
cd Data/Yeast
for dir in ERR*; do
    fn=`awk -v id=$dir '$1 == id { printf("%s-%s-%s.fastq.gz", $1, $3, $4) }' \
	../../ERP004763_sample_mapping.tsv`
    echo $fn
    cd Raw-renamed
    ln -sf ../$dir/$dir.fastq.gz $fn
    cd ..
done
