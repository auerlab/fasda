#!/bin/sh -e

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

gff_filename=Results/04-reference/$(Reference/gff-filename.sh)
ab_dir=Results/10-fasda-hisat2

for bam in Results/09-hisat2-align/*.bam; do
    base=$(basename $bam)
    ab=$ab_dir/${base%.bam}-abundance.tsv
    if [ ! -e $ab ]; then
	fasda abundance --stringtie --output-dir Results/10-fasda-hisat2 \
	    75 $gff_filename $bam
    else
	printf "$ab already exists.\n"
    fi
done

transcripts='ENSMUST00000000127 ENSMUST00000000543 ENSMUST00000019143'

for transcript in $transcripts; do
    printf "\nCondition 1 abundances:\n"
    awk -v t=$transcript '$1 == t { print $1, $4 }' Results/10-fasda-hisat2/cond1-rep*
done
pause

# Now compute fold-changes and P-values for various sets of 3 replicates
# out of the 6 available
for first in $(seq 1 4); do
    for second in $(seq $(($first + 1)) 5); do
	for third in $(seq $(($second + 1)) 6); do
	    for cond in 1 2; do
		if [ ! -e $ab_dir/norm-cond$cond-$first-$second-$third.tsv ]; then
		    printf "Normalizing cond1 replicates $first $second $third...\n"
		    fasda normalize --output \
			$ab_dir/norm-cond$cond-$first-$second-$third.tsv \
			$ab_dir/cond$cond-rep$first-abundance.tsv \
			$ab_dir/cond$cond-rep$second-abundance.tsv \
			$ab_dir/cond$cond-rep$third-abundance.tsv
		else
		    printf "Using existing $ab_dir/norm-$first-$second-$third.tsv.\n"
		fi
	    done
	    
	    #printf "Computing fold-changes for $first-$second-$third...\n"
	    #fasda fold-change --output $ab_dir/fc-$first-$second-$third.txt \
	    #    $ab_dir/cond1-rep$first-abundance.tsv \
	    #    $ab_dir/cond2-rep$first-abundance.tsv
	done
    done
done

for transcript in $transcripts; do
    printf "===\n"
    # awk -v t=$transcript '$1 == t { print FILENAME, $1, $4 }' $ab_dir/norm-cond1-*.tsv
    fgrep $transcript $ab_dir/norm-cond1*.tsv
done


