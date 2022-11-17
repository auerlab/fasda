#!/bin/sh -e

# Document software versions used for publication
uname -a
fastq-trim --version
pwd

# So gzip can keep up with fastq-trim
export GZIP=-1

trimmed_dir=Data/02-trim
mkdir -p $trimmed_dir
for file in Data/Raw-renamed/*.fastq.gz; do
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
