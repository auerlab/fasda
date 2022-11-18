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
    time fasda normalize --output $condition-all-norm-$tr.tsv \
	$kallisto_dir/$condition-*/abundance.tsv
done

printf "Computing fold-change...\n"
time fasda fold-change --output WT-SNF2-FC-MW-$tr.txt \
    WT-all-norm-$tr.tsv SNF2-all-norm-$tr.tsv
head WT-SNF2-FC-MW-$tr.txt

##########################################################################
#   3 to 12 replicates, [near-]exact P-values
##########################################################################

if [ $tr -lt 12 ]; then
    max_ne=$tr
else
    max_ne=12
fi

for replicates in $(seq 3 $max_ne); do
    for condition in WT SNF2; do
	if [ ! -e $condition-all-norm-$replicates.tsv ]; then
	    printf "Normalizing $condition: $replicates replicates\n"
	    time fasda normalize --output \
		$condition-all-norm-$replicates.tsv \
		$kallisto_dir/$condition-[1-$replicates]/abundance.tsv
	fi
    done
    
    if [ ! -e WT-SNF2-FC-NE-$replicates.txt ]; then
	if [ ! -e WT-SNF2-FC-NE-$replicates.txt ]; then
	    printf "Computing fold-change for $replicates replicates...\n"
	    time fasda fold-change --near-exact \
		--output WT-SNF2-FC-NE-$replicates.txt \
		WT-all-norm-$replicates.tsv SNF2-all-norm-$replicates.tsv
	fi
    fi
done

##########################################################################
#   8 to all replicates, Mann-Whitney P-values
##########################################################################

for replicates in $(seq 8 $max_ne) `seq $(($max_ne + 1)) 5 $tr`; do
    for condition in WT SNF2; do
	if [ ! -e $condition-all-norm-$replicates.tsv ]; then
	    printf "Normalizing $condition: $replicates replicates\n"
	    time fasda normalize --output \
		$condition-all-norm-$replicates.tsv \
		$kallisto_dir/$condition-[1-$replicates]/abundance.tsv
	fi
    done
    
    if [ ! -e WT-SNF2-FC-MW-$replicates.txt ]; then
	if [ ! -e WT-SNF2-FC-MW-$replicates.txt ]; then
	    printf "Computing fold-change for $replicates replicates...\n"
	    time fasda fold-change \
		--output WT-SNF2-FC-MW-$replicates.txt \
		WT-all-norm-$replicates.tsv SNF2-all-norm-$replicates.tsv
	fi
    fi
done
    
head WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt | more
printf "%-25s %10s %10s\n" "File" "Features" "P < 0.05"
for file in WT-SNF2-FC-NE-*.txt WT-SNF2-FC-MW-*.txt; do
    printf "%-25s %10s %10s\n" $file: $(cat $file | wc -l) $(awk '$5 < 0.05' $file | wc -l)
done
