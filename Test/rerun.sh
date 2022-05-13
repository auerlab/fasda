#!/bin/sh -e

##########################################################################
#   Description:
#       Remove output files and logs from a previous run and resubmit
#       
#   History:
#   Date        Name        Modification
#   2021-11-27  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 script.sh\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
script=$1

base=${script%.s*h}
printf "Remove results from Data/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Data/$base/*
fi

printf "Remove logs from Logs/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -rf Logs/$base/*
fi

if [ ${script##*.} = sbatch ]; then
    sbatch $script
else
    ./$script
fi
