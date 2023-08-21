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
    fi
done

# Transcripts with modest to high coverage and likely significant
transcripts='ENSMUST00000000090 ENSMUST00000000109 ENSMUST00000029451'

if false; then
for transcript in $transcripts; do
    printf "\nCondition 1 abundances:\n"
    awk -v t=$transcript '$1 == t { print $1, $4 }' Results/10-fasda-hisat2/cond1-rep*
done
pause
fi

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
		fi
	    done
	    
	    if [ ! -e $ab_dir/fc-$first-$second-$third.txt ]; then
		printf "Computing fold-changes for $first-$second-$third...\n"
		fasda fold-change --output $ab_dir/fc-$first-$second-$third.txt \
		    $ab_dir/norm-cond1-$first-$second-$third.tsv \
		    $ab_dir/norm-cond2-$first-$second-$third.tsv
	    fi
	done
    done
done

# Compute DE for all 6 replicates
for cond in 1 2; do
    if [ ! -e $ab_dir/norm-cond$cond-all.tsv ]; then
	printf "Normalizing condition $cond for all replicates...\n"
	fasda normalize --output \
	    $ab_dir/norm-cond$cond-all.tsv $ab_dir/cond$cond-rep*-abundance.tsv
    fi
done
if [ ! -e $ab_dir/fc-all.txt ]; then
    printf "Computing fold-changes for all replicates...\n"
    fasda fold-change --output $ab_dir/fc-all.txt \
	$ab_dir/norm-cond1-all.tsv $ab_dir/norm-cond2-all.tsv
fi

cd $ab_dir
for transcript in $transcripts; do
    printf "\n=== $transcript ===\n"
    awk -v t=$transcript \
	'$1 == t { printf("%-20s %10.1f %10.1f %10.1f\n", FILENAME, $2, $3, $4); }' \
	norm-cond1-1-2-3.tsv norm-cond1-4-5-6.tsv
done

for transcript in $transcripts; do
    printf "\n=== $transcript ===\n"
    printf "%-15s %6s %6s %4s %4s %4s %5s %4s\n" \
	"Samples" "MNC1" "MNC2" "SD1" "SD2" "Agr" "FC" "P"
    awk -v t=$transcript \
	'$1 == t { printf("%-15s %6.1f %6.1f %4.1f %4.1f %4d %5.1f %4.2f\n",
	    FILENAME, $2, $3, $4, $5, $6, $7, $8); }' \
	*.txt | sort -k 8 -n
done

