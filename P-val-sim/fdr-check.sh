#!/bin/sh -e

./pval 200 200 .3 5 10000 \
    | awk '$5 == "P-value" && $11 <= 0.01 { print $11 }' | wc -l

