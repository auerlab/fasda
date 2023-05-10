#!/bin/sh -e

##########################################################################
#   Description:
#       Compute fold-changes and P-values for N randomly chosen samples
#       from the 48 available for the yeast data.  This demonstrates how
#       variable fold-changes and P-values can be due to biological
#       and technical variation in the coverages.
#       
#   History:
#   Date        Name        Modification
#   2023-05-10  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 sample-count\n"
    exit 1
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2023-05-10  Jason Bacon Begin
##########################################################################

random()
{
    if [ $# != 1 ]; then
	printf "Usage: random sample-count\n"
	exit 1
    fi
    modulus=$1
    
    taken=true
    while [ $taken = true ]; do
	random=$(head /dev/urandom | cksum | cut -d ' ' -f 1)
	random=$(($random % $modulus))
	if [ ! -e taken.$random ]; then
	    taken=false
	    touch taken.$random
	fi
    done
    printf "$random\n"
}

##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
sample_count=$1

kallisto_dir=../06-kallisto-quant
cd Results/07-fasda-kallisto

# Compare 10 random selections of samples
for trial in $(seq 1 20); do
    printf "\n=== Trial $trial ===\n"
    
    # Choose $sample_count biological samples at random from the 48
    # yeast samples
    rm -f taken.*
    samples=""
    for n in $(seq 1 $sample_count); do
	sample=$(random 48)         # 0-based
	sample=$(($sample + 1))     # 1-based
	wt=$kallisto_dir/WT-$sample/abundance.tsv
	snf2=$kallisto_dir/SNF2-$sample/abundance.tsv
	
	if [ ! -e $wt ] || [ ! -e $snf2 ]; then
	    printf "Error: Missing abundances for sample $sample.\n"
	    printf "Run kallisto on all 48 samples first.\n"
	    exit 1
	fi
	samples="$samples $sample"
    done
    rm -f taken.*
    
    # Compute fold-changes and P-values for this set of samples
    for condition in WT SNF2; do
	r0=$(printf '%02d' $trial)
	if [ ! -e $condition-all-norm-$r0.tsv ]; then
	    printf "Normalizing $condition: $samples replicates\n"
	    files=""
	    for r in $samples; do
		file=$kallisto_dir/$condition-$r/abundance.tsv
		if [ ! -e $file ]; then
		    printf "Missing $file.\n"
		    exit 1
		fi
		files="$files $file"
	    done
	    printf "%s\n" $files
	    time fasda normalize --output $condition-all-norm-$r0.tsv $files
	fi
    done
    
    if [ ! -e WT-SNF2-FC-MW-$r0.txt ]; then
	printf "Computing fold-change for $replicates replicates...\n"
	time fasda fold-change \
	    --output WT-SNF2-FC-MW-$r0.txt \
	    WT-all-norm-$r0.tsv SNF2-all-norm-$r0.tsv
    fi
done

transcripts=$(fgrep -v Feature WT-SNF2-FC-MW-01.txt | head | awk '{ print $1 }')
for transcript in $transcripts; do
    printf "\n$transcript\n"
    awk -v transcript=$transcript '$1 == transcript { print $7, $8 }' *.txt
done
