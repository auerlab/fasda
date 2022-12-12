#!/bin/sh -e

cd ..
make
cd Test

if [ ! -e Data/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 ]; then
    Reference/fetch-gff.sh
fi

for condition in WT SNF2; do
    for replicate in 1 2 3; do
	dir=Data/08-hisat2-align
	file=$condition-$replicate.bam
	ab=$dir/${file%.bam}-abundance.tsv
	if [ ! -e $ab ]; then
	    printf "Computing depth for $condition replicate $replicate...\n"
	    ../abundance \
		Data/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
		$dir/$file
	    head $ab
	fi
    done

    ../normalize --output hisat-$condition-normalized.tsv \
	$dir/$condition-[123]-abundance.tsv
done

../fold-change hisat-WT-normalized.tsv hisat-SNF2-normalized.tsv | more
