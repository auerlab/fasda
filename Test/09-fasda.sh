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

cd Data/09-fasda

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH
which fasda

kallisto_dir=../06-kallisto-quant

##########################################################################
#   All replicates, should use Mann-Whitney
##########################################################################

for condition in WT SNF2; do
    r=$(ls $kallisto_dir/$condition-*/abundance.tsv | wc -l)
    tr=$(echo $r)
    printf "Normalizing $condition: $tr replicates\n"
    time fasda normalize --output $condition-all-norm.tsv \
	$kallisto_dir/$condition-*/abundance.tsv
done

printf "Computing fold-change...\n"
time fasda fold-change --output WT-SNF2-FC-MW-$tr.txt \
    WT-all-norm.tsv SNF2-all-norm.tsv
pause
head WT-SNF2-FC-MW-$tr.txt

printf "Compute near-exact P-values for $tr replicates? y/[n] "
read ne
if [ 0$ne = 0y ]; then
    time fasda fold-change --near-exact --output WT-SNF2-FC-NE-$tr.txt \
	WT-all-norm.tsv SNF2-all-norm.tsv
fi
if [ -e WT-SNF2-FC-NE-$tr.txt ]; then
    head WT-SNF2-FC-MW-$tr.txt WT-SNF2-FC-NE-$tr.txt
    pause
fi

##########################################################################
#   3 - 7 replicates, should use exact P-values
##########################################################################

for replicates in $(seq 3 7); do
    for condition in WT SNF2; do
	printf "Normalizing $condition: $replicates replicates\n"
	time fasda normalize --output \
	    $condition-all-norm-$replicates.tsv \
	    $kallisto_dir/$condition-[1-$replicates]/abundance.tsv
    done
    
    printf "Computing fold-change for $replicates replicates...\n"
    time fasda fold-change --output WT-SNF2-FC-NE-$replicates.txt \
	WT-all-norm-$replicates.tsv SNF2-all-norm-$replicates.tsv
    pause
    
    head WT-SNF2-FC-MW-$tr.txt WT-SNF2-FC-NE-$tr.txt WT-SNF2-FC-NE-*.txt | more
    pause
done
