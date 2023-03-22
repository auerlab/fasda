#!/bin/sh -e

printf "Are you sure you want to remove all results? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Results/0[2-9]* Results/1[0-9]*
fi
ls Results
