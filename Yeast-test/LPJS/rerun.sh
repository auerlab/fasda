#!/bin/sh -e

##########################################################################
#   Description:
#       Remove output files and logs from a previous run and resubmit
##########################################################################

usage()
{
    printf "Usage: $0 script\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
script=$1

base=${script%.*}
printf "Remove results from Results/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Results/$base/*
fi

log_dir=LPJS-logs/$base
printf "Remove logs from $log_dir? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf $log_dir/*
fi

if [ ${script##*.} = sbatch ]; then
    sbatch $script
elif [ ${script##*.} = lpjs ]; then
    lpjs submit $script
else
    ./$script
fi
