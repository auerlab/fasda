#!/bin/sh -e

for p_val in 0.05 0.01; do
    case $p_val in
    0.05)
	runs=1000
	;;
    0.01)
	runs=10000
	;;
    esac
    for replicates in $(seq 3 7); do
	sum=0
	for i in $(seq 1 10); do
	    count=`./pval 100 100 .2 $replicates $runs \
		| awk '$5 == "P-value" && $11 <= '$p_val' { print $11 }' \
		| wc -l`
	    sum=$((sum + count))
	done
	printf "P-value: $p_val  Replicates: $replicates  "
	printf "Average of $i x $runs trials: $((sum / i))\n"
    done
done
