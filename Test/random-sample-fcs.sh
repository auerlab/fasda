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

# Compare 10 random selections of samples
for trial in $(seq 1 10); do

    # Choose $sample_count biological samples at random from the 48
    # yeast samples
    rm -f taken.*
    for n in $(seq 1 $sample_count); do
	sample=$(random 48)         # 0-based
	sample=$(($sample + 1))     # 1-based
	wt=Results/06-kallisto-quant/WT-$sample/abundance.tsv
	snf2=Results/06-kallisto-quant/SNF2-$sample/abundance.tsv
	echo $wt $snf2
    done
    rm -f taken.*
    
    # Compute fold-changes and P-values for this set of samples
    
done
