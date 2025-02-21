#!/bin/sh -e

##########################################################################
#   Description:
#       Run fasda normalize and fold-change on star abundances
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

cd Results/13-fasda-star

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH

uname -a
fasda --version
pwd

reference_dir=../04-reference
star_dir=../12-star-align
log_dir=../../Logs/13-fasda-star

##########################################################################
#   Compute abundances
##########################################################################

for condition in WT SNF2; do
    for r in $(seq 1 $tr); do
	file=$condition-$r/Aligned.sortedByCoord.out.bam
	ab=$star_dir/${file%.bam}-abundance.tsv
	printf "Computing abundances for $condition replicate $r...\n"
	set -x
	time fasda abundance 50 \
	    $reference_dir/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
	    $star_dir/$file
	set +x
	column -t $ab | head
	wc $ab
    done
done
exit

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
printf "Hyperthreads = $threads  Jobs = $jobs\n"

seq 3 $max_ne | xargs -n 1 -P $jobs \
    ../../fasda-star-ne.sh $star_dir $log_dir

exit

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
		    files="$files $star_dir/$condition-$r/abundance.tsv"
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

if [ $tr -ge 8 ]; then
    head WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt | more
fi

printf "\n%-25s %10s %10s\n" "File" "Features" "P < 0.05"
for file in WT-SNF2-FC-NE-*.txt; do
    printf "%-25s %10s %10s\n" $file: \
	$(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
done | more

if [ $tr -ge 8 ] && [ -n "$(ls WT-SNF2-FC-MW-*.txt)" ]; then
    for file in WT-SNF2-FC-MW-*.txt; do
	printf "%-25s %10s %10s\n" $file: \
	    $(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
    done | more
fi
