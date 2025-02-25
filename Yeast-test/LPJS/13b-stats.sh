#!/bin/sh -e

if [ $# = 1 ]; then
    results_dir=$1
else
    results_dir=Results
fi

printf "Total features:\n"
for file in $results_dir/13-fasda-kallisto/*.txt; do
    fgrep -v P-val $file > $file-no-header
    wc -l $file-no-header
    rm -f $file-no-header
done
total=$(wc -l $results_dir/13-fasda-kallisto/*.txt | awk '{ print $1 }')

printf "\nFeatures with P-values < 0.05:\n"
for file in $results_dir/13-fasda-kallisto/*.txt; do
    sig=$(awk '$8 < 0.05' $file | wc -l | awk '{ print $1 }')
    percent=$(printf "scale=2\n100 * $sig / $total\n" | bc)
    printf "$(awk '$8 < 0.05' $file | wc -l) ($percent%%) $file\n"
done
