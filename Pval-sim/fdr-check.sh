#!/bin/sh -e

make
for p_val in 0.05 0.01; do
    case $p_val in
    0.10|0.05)
	trials=1000
	;;
    0.01)
	trials=10000
	;;
    esac
    rm -f tmp*.txt
    for replicates in $(seq 3 10); do
	printf "P-value: $p_val  Replicates: $replicates  "
	printf "Expected false positives = %s\n" \
	    $(printf "$p_val * $trials\nquit\n" | bc -l)
	
	# Run expensive P-value computations in parallel
	seq 1 4 | parallel ./pval-sim 100 100 .5 $replicates $trials '>' tmp{}.txt
	
	for run in $(seq 1 4); do
	    awk '$1 == "Exact" && $4 <= '$p_val' { print $4 }' tmp$run.txt \
		| wc -l > exact-output.$run
	    awk '$1 == "Exact" && $7 <= '$p_val' { print $7 }' tmp$run.txt \
		| wc -l > mw-output.$run
	done
	
	exact_counts=$(cat exact-output.*)
	exact_sum=0
	for v in $exact_counts; do
	    exact_sum=$((exact_sum + v))
	done
	printf "Exact P-val counts for each run of $trials trials: "
	printf "%s " $exact_counts
	printf "\nExact P-val average of $run x $trials trials: $(($exact_sum / $run))\n"
	
	if [ $replicates -ge 8 ]; then
	    mw_counts=$(cat mw-output.*)
	    mw_sum=0
	    for v in $mw_counts; do
		mw_sum=$((mw_sum + v))
	    done
	    printf "Mann-Whitney counts for each run of $trials trials: "
	    printf "%s " $mw_counts
	    printf "\nMann-Whitney average of $run x $trials trials: $(($mw_sum / $run))\n"
	fi
	
	printf "\n"
    done
    rm -f tmp*.txt exact-output.* mw-output.*
done
