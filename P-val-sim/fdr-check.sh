#!/bin/sh -e

make
for p_val in 0.05 0.01; do
    case $p_val in
    0.10|0.05)
	runs=1000
	;;
    0.01)
	runs=10000
	;;
    esac
    rm -f tmp*.txt
    for replicates in $(seq 3 10); do
	# Run expensive P-value computations in parallel
	seq 1 4 | parallel ./pval 100 100 .5 $replicates $runs '>' tmp{}.txt
	
	for i in $(seq 1 4); do
	    awk '$1 == "Exact" && $4 <= '$p_val' { print $4 }' tmp$i.txt \
		| wc -l > exact-output.$i
	    awk '$1 == "Exact" && $7 <= '$p_val' { print $7 }' tmp$i.txt \
		| wc -l > mw-output.$i
	done
	
	exact_counts=$(cat exact-output.*)
	exact_sum=0
	for v in $exact_counts; do
	    exact_sum=$((exact_sum + v))
	done

	mw_counts=$(cat mw-output.*)
	mw_sum=0
	for v in $mw_counts; do
	    mw_sum=$((mw_sum + v))
	done
	
	printf "Expected false positives = %s\n" \
	    $(printf "$p_val * $runs\nquit\n" | bc -l)
	printf "P-value: $p_val  Replicates: $replicates  "
	printf "Exact average of $i x $runs trials: $((exact_sum / i))\n"
	printf "Exact counts for each run: "
	printf "%s " $exact_counts
	printf "\nMann-whitney average of $i x $runs trials: $((mw_sum / i))\n"
	printf "Mann-Whitney counts for each run: "
	printf "%s " $mw_counts
	printf "\n\n"
    done
    rm -f tmp*.txt exact-output.* mw-output.*
done
