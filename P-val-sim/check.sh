#!/bin/sh -e

run_1000()
{
    for c in $(seq 1 1000); do
	./pval 200 200 .3 3 | fgrep 'FC count' \
	    | awk '$14 <= 0.05 { print $14 }'
    done
}

run_1000 | wc -l
