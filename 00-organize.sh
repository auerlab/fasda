#!/bin/sh -e

mkdir -p Data Logs
scripts=$(ls 0[2-9]-*) # [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Data/$stage Logs/$stage
done

mkdir -p Data/Raw-renamed
cd Data/Raw-renamed
for dir in ../../Yeast/ERR*; do
    base=$(basename $dir)
    test -e $dir/$base.fastq.gz
    fn=`awk -v id=$base '$1 == id { printf("%s-%s.fastq.gz", $3, $4) }' \
	../../ERP004763_sample_mapping.tsv`
    printf "$dir/$base.fastq.gz -> $fn\n"
    ln -sf $dir/$base.fastq.gz $fn
done
