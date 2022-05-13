#!/bin/sh -e

# 2020-10-21 time ./proximal-and-distal-peaks.sh
# nawk (BSD/Mac): 35 seconds
# mawk: 19 seconds
# gawk: 20 seconds
for awk in mawk gawk awk; do
    if which $awk > /dev/null 2>&1; then
	echo $awk
	exit
    fi
done
    
