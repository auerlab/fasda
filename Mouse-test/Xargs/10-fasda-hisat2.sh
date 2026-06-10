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

cd ..
./cave-man-install.sh
cd Mouse-test
PATH=../../local/bin:${PATH}
export PATH

# FIXME: stringtie requires a GTF for gene-level abundances, while
# our tools use gff3
gff_filename=Results/04-reference/$(Reference/gtf-filename.sh)
ab_dir=Results/10-fasda-hisat2
mkdir -p $ab_dir

for bam in Results/09-hisat2-align/*.bam; do
    base=$(basename $bam)
    ab=$ab_dir/${base%.bam}-abundance.tsv
    if [ ! -e $ab ]; then
	fasda abundance --stringtie --output-dir Results/10-fasda-hisat2 \
	    75 $gff_filename $bam
    fi
done

# Now compute fold-changes and P-values for various sets of 3 replicates
# out of the 6 available
for first in $(seq 1 4); do
    for second in $(seq $(($first + 1)) 5); do
	for third in $(seq $(($second + 1)) 6); do
	    norm_file=$ab_dir/norm-all-$first-$second-$third.tsv
	    printf "Generating $norm_file...\n"
	    set -x
	    fasda normalize --output $norm_file \
		$ab_dir/*-rep$first-abundance.tsv \
		$ab_dir/*-rep$second-abundance.tsv \
		$ab_dir/*-rep$third-abundance.tsv
	    set +x
	    head $norm_file
	    
	    cut -f 1-4 $norm_file > $ab_dir/norm-cond1-$first-$second-$third.tsv
	    cut -f 1,5-7 $norm_file > $ab_dir/norm-cond2-$first-$second-$third.tsv
	    printf "Computing fold-changes for $first-$second-$third...\n"
	    fasda fold-change --output $ab_dir/fc-$first-$second-$third.txt \
		$ab_dir/norm-cond1-$first-$second-$third.tsv \
		$ab_dir/norm-cond2-$first-$second-$third.tsv
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
