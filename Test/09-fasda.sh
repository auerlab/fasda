#!/bin/sh -e

##########################################################################
#   Description:
#       Run fasda normalize and fold-change on kallisto abundances
#       
#   History:
#   Date        Name        Modification
#   2022-11-19  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 max-replicates\n"
    exit 1
}


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
tr=$1

cd Data/09-fasda

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH
which fasda

kallisto_dir=../06-kallisto-quant
log_dir=../../Logs/09-fasda

##########################################################################
#   3 to 12 replicates, [near-]exact P-values
#   Run in parallel since some can take a few minutes
##########################################################################

if [ $tr -lt 12 ]; then
    max_ne=$tr
else
    max_ne=12
fi

if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi
jobs=$threads
if [ $jobs = 0 ]; then
    jobs=1
fi
printf "Hyperthreads = $threads  Jobs = $jobs\n"

seq 3 $max_ne | xargs -n 1 -P $jobs \
    ../../fasda-ne.sh $kallisto_dir $log_dir

##########################################################################
#   8 to all replicates, Mann-Whitney P-values
##########################################################################

if [ $tr -ge 8 ]; then
    for replicates in $(seq 8 $tr); do
	# FIXME: Factor out to fasda-mw.sh?
	r0=$(printf '%02d' $replicates)
	for condition in WT SNF2; do
	    if [ ! -e $condition-all-norm-$r0.tsv ]; then
		printf "Normalizing $condition: $replicates replicates\n"
		files=""
		for r in $(seq 1 $replicates); do
		    files="$files $kallisto_dir/$condition-$r/abundance.tsv"
		done
		# printf "%s\n" $files
		time fasda normalize --output \
		    $condition-all-norm-$r0.tsv $files \
		    > $log_dir/normalize-$condition-$r0-MW.out \
		    2> $log_dir/normalize-$condition-$r0-MW.err
	    fi
	done
	
	if [ ! -e WT-SNF2-FC-MW-$r0.txt ]; then
	    printf "Computing fold-change for $replicates replicates...\n"
	    time fasda fold-change \
		--output WT-SNF2-FC-MW-$r0.txt \
		WT-all-norm-$r0.tsv SNF2-all-norm-$r0.tsv \
		> $log_dir/fc-$condition-$r0-MW.out \
		2> $log_dir/fc-$condition-$r0-MW.err
	fi
    done
fi

head WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt | more
printf "\n%-25s %10s %10s\n" "File" "Features" "P < 0.05"
for file in WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt; do
    printf "%-25s %10s %10s\n" $file: \
	$(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
done | more
