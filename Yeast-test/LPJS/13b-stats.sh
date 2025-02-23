#!/bin/sh -e

printf "Total features:\n"
wc -l Results/13-fasda-kallisto/*.txt
total=$(wc -l Results/13-fasda-kallisto/*.txt | awk '{ print $1 }')

printf "\nFeatures with P-values < 0.05:\n"
for file in Results/13-fasda-kallisto/*.txt; do
    sig=$(awk '$8 < 0.05' $file | wc -l | awk '{ print $1 }')
    percent=$(printf "scale=2\n100 * $sig / $total\n" | bc)
    printf "$(awk '$8 < 0.05' $file | wc -l) ($percent%%) $file\n"
done
