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
cd Test

./00-organize.sh

if [ ! -e Data/Raw-renamed ]; then
    printf "How many samples would you like to download? [10] "
    read count
    if [ -z "$count" ]; then
	count=10
    fi
    ./01-fetch.sh $count
fi
raw_count=$(ls Data/Raw-renamed/*.gz | wc -l)

trimmed_count=$(ls Data/02-trim/*.gz | wc -l)
if [ $trimmed_count -ne $raw_count ]; then
    ./02-trim.sh
else
    printf "02-trim.sh already done.\n"
fi

qc_count=$(ls Data/03-qc/02-trim/*.zip | wc -l)
if [ $qc_count -ne $raw_count ]; then
    ./03-qc.sh
else
    printf "03-qc.sh already done.\n"
fi

if [ ! -e Data/04-reference/all-but-xy.genome.fa.fai ]; then
    ./04-reference.sh
else
    printf "04-reference.sh already done.\n"
fi

if [ ! -e Data/05-kallisto-index/all-but-xy.index ]; then
    ./05-kallisto-index.sh
else
    printf "05-kallisto-index.sh already done.\n"
fi

quant_count=$(ls -d Data/06-kallisto-quant/* | wc -l)
if [ $quant_count -ne $raw_count ]; then
    ./06-kallisto-quant.sh
else
    printf "06-kallisto-quant.sh already done.\n"
fi

dir=Data/06-kallisto-quant
PATH=../../local/bin:$PATH
export PATH
which fasda

##########################################################################
#   All replicates, should use Mann-Whitney
##########################################################################

for condition in WT SNF2; do
    printf "Normalizing $condition...\n"
    time fasda normalize --output $condition-all-norm.tsv \
	$dir/$condition-*/abundance.tsv
done

printf "Computing fold-change...\n"
time fasda fold-change --output WT-SNF2-FC-MW.txt \
    WT-all-norm.tsv SNF2-all-norm.tsv
pause
more WT-SNF2-FC-MW.txt

printf "Compute near-exact P-values for all replicates? y/[n] "
read ne
if [ 0$ne = 0y ]; then
    which fasda
    time fasda fold-change --near-exact --output WT-SNF2-FC-NE.txt \
	WT-all-norm.tsv SNF2-all-norm.tsv
    pause
    more WT-SNF2-FC-NE.txt
fi
if [ -e WT-SNF2-FC-NE.txt ]; then
    head WT-SNF2-FC-MW.txt WT-SNF2-FC-NE.txt
    pause
fi

##########################################################################
#   3 - 7 replicates, should use exact P-values
##########################################################################

for replicates in $(seq 3 7); do
    for condition in WT SNF2; do
	printf "Normalizing $condition...\n"
	time fasda normalize --output \
	    $condition-all-norm-$replicates.tsv \
	    $dir/$condition-[1-$replicates]/abundance.tsv
    done
    
    printf "Computing fold-change for $replicates replicates...\n"
    time fasda fold-change --output WT-SNF2-FC-NE-$replicates.txt \
	WT-all-norm-$replicates.tsv SNF2-all-norm-$replicates.tsv
    pause
    more WT-SNF2-FC-NE-$replicates.txt
    
    head WT-SNF2-FC-MW.txt WT-SNF2-FC-NE.txt WT-SNF2-FC-NE-*.txt | more
    pause
done
