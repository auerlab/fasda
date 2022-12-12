#!/bin/sh -e

##########################################################################
#   Description:
#       Compute fold-change from hisat2 data
#       
#   History:
#   Date        Name        Modification
#   2022-12-11  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 max-replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
max_replicates=$1

cd ..
make
cd Test

if [ ! -e Data/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 ]; then
    Reference/fetch-gff.sh
fi

hisat_dir=Data/09-hisat2-align
fasda_dir=Data/10-fasda-hisat2
for condition in WT SNF2; do
    for r in $(seq 1 $max_replicates); do
	file=$condition-$r.bam
	ab=$hisat_dir/${file%.bam}-abundance.tsv
	if [ ! -e $ab ]; then
	    printf "Computing depth for $condition replicate $r...\n"
	    ../abundance \
		Data/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
		$hisat_dir/$file
	    head $ab
	fi
    done
done

for replicates in $(seq 2 $max_replicates); do
    for condition in WT SNF2; do
	rep_files=""
	for r in $(seq 1 $replicates); do
	    rep_files="$rep_files $hisat_dir/$condition-$r-abundance.tsv"
	done
	printf "Normalizing $rep_files...\n"
	../normalize --output $fasda_dir/$condition-all-norm-$replicates.tsv \
	    $rep_files
    done
    
    printf "Computing fold-change for $replicates replicates...\n"
    fc_out=$fasda_dir/WT-SNF2-FC-MW-$replicates.txt
    ../fold-change \
	--output $fc_out \
	$fasda_dir/WT-all-norm-$replicates.tsv \
	$fasda_dir/SNF2-all-norm-$replicates.tsv
    
    for feature in YPL071C_mRNA YLL050C_mRNA YMR172W_mRNA YOR185C_mRNA; do
	grep $feature $fc_out
    done
done
