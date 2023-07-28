#!/bin/sh -e

##########################################################################
#   Description:
#       Run fasda normalize and fold-change on hisat abundances
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

uname -a
fasda --version
pwd

if [ ! -e Results/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 ]; then
    Reference/fetch-gff.sh
fi

cd Results/10-fasda-hisat

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH

hisat_dir=../09-hisat-align
reference_dir=../04-reference
log_dir=../../Logs/10-fasda-hisat

##########################################################################
#   3 to 12 replicates, [near-]exact P-values
#   Run in parallel since some can take a few minutes
##########################################################################

if [ $tr -lt 12 ]; then
    max_ne=$tr
else
    max_ne=12
fi

##########################################################################
#   Compute abundances
##########################################################################

for condition in WT SNF2; do
    for r in $(seq 1 $tr); do
	file=$condition-$r.bam
	ab=$hisat_dir/${file%.bam}-abundance.tsv
	    printf "Computing abundances for $condition replicate $r...\n"
	    time fasda abundance 50 \
		$reference_dir/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
		$hisat_dir/$file
	    
	    # FIXME: Hack for testing, move this to abundance.c
	    # stringtie is screwy, sorting each gtf output differently
	    #fgrep eff_length $ab > temp.tsv
	    #fgrep -v eff_length $ab | sort >> temp.tsv
	    #mv -f temp.tsv $ab
	    
	    column -t $ab | head
	    wc $ab
    done
done

if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi
jobs=$threads
printf "Hyperthreads = $threads  Jobs = $jobs\n"

seq 3 $max_ne | xargs -n 1 -P $jobs \
    ../../fasda-hisat-ne.sh $hisat_dir $log_dir

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
		    files="$files $hisat_dir/$condition-$r-abundance.tsv"
		done
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
for file in WT-SNF2-FC-NE-*.txt; do
    printf "%-25s %10s %10s\n" $file: \
	$(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
done | more

if [ -n "$(ls WT-SNF2-FC-MW-*.txt)" ]; then
    for file in WT-SNF2-FC-MW-*.txt; do
	printf "%-25s %10s %10s\n" $file: \
	    $(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
    done | more
fi

printf "\nHisat:\n"
for feature in YPL071C_mRNA YLL050C_mRNA YMR172W_mRNA YOR185C_mRNA; do
    grep -h $feature WT-SNF2-FC-NE-03.txt
done
printf "\nKallisto:\n"
for feature in YPL071C_mRNA YLL050C_mRNA YMR172W_mRNA YOR185C_mRNA; do
    grep -h $feature ../07-fasda-kallisto/WT-SNF2-FC-NE-03.txt
done

