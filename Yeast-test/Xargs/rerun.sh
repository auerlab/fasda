#!/bin/sh -e

##########################################################################
#   Description:
#       Remove output files and logs from a previous run and resubmit
##########################################################################

usage()
{
    printf "Usage: $0 script.sh [args]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi
script=$1
shift

base=${script%.s*h}
printf "Remove results from Results/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Results/$base/*
fi

printf "Remove logs from Logs/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Logs/$base/*
fi

if [ ${script##*.} = sbatch ]; then
    sbatch $script "$@"
else
    ./$script "$@"
fi
