#!/bin/sh -e

trimmed_dir=Data/01-trim
mkdir -p $trimmed_dir
for file in Data/Yeast/Raw-renamed/*.fastq.gz; do
    base=`basename $file`
    trimmed=$trimmed_dir/$base
    if [ -e $trimmed ]; then
	printf "$trimmed already exists.\n"
    else
	# Adapter discovered by fastq-scum
	fastq-trim --3p-adapter1 AGATCGGAAGAG --polya-min-length 3 \
	    $file $trimmed
    fi
done
