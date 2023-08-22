#!/bin/sh -e

##########################################################################
#   Description:
#       Compute one FASDA [near-]exact set of FCs
#       
#   History:
#   Date        Name        Modification
#   2022-11-27  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 star-dir replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
star_dir=$1
log_dir=$2
replicates=$3

r0=$(printf '%02d' $replicates)
for condition in WT SNF2; do
    if [ ! -e $condition-all-norm-$r0.tsv ]; then
	printf "Normalizing $condition: $replicates replicates\n"
	files=""
	for r in $(seq 1 $replicates); do
	    files="$files $star_dir/$condition-$r/*-abundance.tsv"
	done
	# printf "%s\n" $files
	time fasda normalize --output $condition-all-norm-$r0.tsv $files \
	    > $log_dir/normalize-$condition-$r0-NE.out \
	    2> $log_dir/normalize-$condition-$r0-NE.err
    fi
done

if [ ! -e WT-SNF2-FC-NE-$r0.txt ]; then
    printf "Computing fold-change for $replicates replicates...\n"
    time fasda fold-change --near-exact \
	--output WT-SNF2-FC-NE-$r0.txt \
	WT-all-norm-$r0.tsv SNF2-all-norm-$r0.tsv \
	> $log_dir/fc-$condition-$r0-NE.out \
	2> $log_dir/fc-$condition-$r0-NE.err
fi

