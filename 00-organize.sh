#!/bin/sh -e

mkdir -p Data Logs
scripts=$(ls 0[2-9]-*) # [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Data/$stage Logs/$stage
done

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
