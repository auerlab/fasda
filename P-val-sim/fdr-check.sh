#!/bin/sh -e

run()
{
    for c in $(seq 1 1000); do
	./pval 200 200 .3 3 \
	    | awk '$5 == "P-value" && $11 <= 0.05 { print $11 }'
    done
}

run | wc -l
