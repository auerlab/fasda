#!/bin/sh -e

##########################################################################
#   Description:
#       Trim a single raw fastq file, storing results in data-dir and
#       redirecting stdout and stderr to log-dir.  This script takes
#       the filename last, so that it can easily be used with xargs
#       to parallelize trimming.
##########################################################################

usage()
{
    printf "Usage: $0 data-dir log-dir file\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
output_dir=$1
log_dir=$2
file=$3

base=`basename $file`
output=$output_dir/$base
if [ -e $output ]; then
    printf "$output already exists.\n"
else
    # Adapter discovered by fastq-scum
    gzcat $file | fastqc -o Results/03-qc stdin:$base
fi
