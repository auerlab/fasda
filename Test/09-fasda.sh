#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       
#   Arguments:
#       
#   Returns:
#
#   Examples:
#
#   Files:
#
#   Environment:
#
#   See also:
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

##########################################################################
#   3 to 12 replicates, [near-]exact P-values
##########################################################################

if [ $tr -lt 12 ]; then
    max_ne=$tr
else
    max_ne=12
fi

for replicates in $(seq 3 $max_ne); do
    r0=$(printf '%02d' $replicates)
    printf "r0 = $r0\n"
    for condition in WT SNF2; do
	if [ ! -e $condition-all-norm-$r0.tsv ]; then
	    printf "Normalizing $condition: $replicates replicates\n"
	    files=""
	    for r in $(seq 1 $replicates); do
		files="$files $kallisto_dir/$condition-$r/abundance.tsv"
	    done
	    printf "%s\n" $files
	    time fasda normalize --output $condition-all-norm-$r0.tsv $files
	fi
    done
    
    if [ ! -e WT-SNF2-FC-NE-$r0.txt ]; then
	printf "Computing fold-change for $replicates replicates...\n"
	time fasda fold-change --near-exact \
	    --output WT-SNF2-FC-NE-$r0.txt \
	    WT-all-norm-$r0.tsv SNF2-all-norm-$r0.tsv
    fi
done

##########################################################################
#   8 to all replicates, Mann-Whitney P-values
##########################################################################

if [ $tr -ge 8 ]; then
    for replicates in $(seq 8 $tr); do
	r0=$(printf '%02d' $replicates)
	printf "r0 = $r0\n"
	for condition in WT SNF2; do
	    if [ ! -e $condition-all-norm-$r0.tsv ]; then
		printf "Normalizing $condition: $replicates replicates\n"
		files=""
		for r in $(seq 1 $replicates); do
		    files="$files $kallisto_dir/$condition-$r/abundance.tsv"
		done
		printf "%s\n" $files
		time fasda normalize --output \
		    $condition-all-norm-$r0.tsv $files
	    fi
	done
	
	if [ ! -e WT-SNF2-FC-MW-$r0.txt ]; then
	    printf "Computing fold-change for $replicates replicates...\n"
	    time fasda fold-change \
		--output WT-SNF2-FC-MW-$r0.txt \
		WT-all-norm-$r0.tsv SNF2-all-norm-$r0.tsv
	fi
    done
fi

head WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt | more
printf "%-25s %10s %10s\n" "File" "Features" "P < 0.05"
for file in WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt; do
    printf "%-25s %10s %10s\n" $file: $(cat $file | wc -l) $(awk '$5 < 0.05' $file | wc -l)
done
