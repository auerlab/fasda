#!/bin/sh -e

cd Results/06-kallisto-quant

# First ten transcripts
genes=$(head -n 11 WT-1/abundance.tsv | tail -10 | cut -f 1)

# All transcripts
# genes=$(cat SNF2-1/abundance.tsv | cut -f 1)

printf "%s\t%6s\t%6s\t%6s\t%7s\n" "Transcript" "Min" "Median" "Max" "Max/Min"
for gene in $genes; do
    awk -v gene=$gene '$1 == gene' WT-*/abundance.tsv | cut -f 5 | sort -g > temp
    # Use pipe to eliminate filename from output
    samples=$(cat temp | wc -l)
    odd=$(($samples % 2))
    median_location=$(($samples / 2))
    # echo $samples $median_location $odd
    min=$(head -1 temp)
    if [ $odd = 1 ]; then
	median=$(head -$median_location temp | tail -1)
    else
	m1=$(head -$median_location temp | tail -1)
	m2=$(head -$(($median_location + 1)) temp | tail -1)
	median=$(echo "scale=2; ($m1 + $m2) / 2" | bc -l)
    fi
    max=$(tail -1 temp)
    if [ $min != 0 ]; then
	spread=$(echo "scale=1; $max / $min" | bc -l)
	printf "%s\t%6.1f\t%6.1f\t%6.1f\t%7.1f\n" $gene $min $median $max $spread
    fi
done
