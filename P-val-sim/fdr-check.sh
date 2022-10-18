#!/bin/sh -e

make
for p_val in 0.10 0.05 0.01; do
    case $p_val in
    0.10|0.05)
	runs=1000
	;;
    0.01)
	runs=10000
	;;
    esac
    for replicates in $(seq 3 7); do
	for i in $(seq 1 4); do
	    ./pval 100 100 .5 $replicates $runs \
		| awk '$5 == "P-value" && $11 <= '$p_val' { print $11 }' \
		| wc -l > output.$i &
	    sum=$((sum + count))
	done
	wait
	counts=$(cat output.*)
	rm -f output.*
	sum=0
	for v in $counts; do
	    sum=$((sum + v))
	done
	printf "P-value: $p_val  Replicates: $replicates  "
	printf "Average of $i x $runs trials: $((sum / i))\n"
	printf "%s " $counts
	printf "\n\n"
    done
done
