#!/bin/sh -e

printf "Are you sure you want to remove all results? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Data/0[2-9]* Data/1[0-9]*
fi
ls Data
