#!/bin/sh -e

##########################################################################
#   Description:
#       Run 1 file through fastqc
#       
#   History:
#   Date        Name        Modification
#   2022-11-28  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 \n"
    exit 1
}

##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
report_dir=$1
log_dir=$2
file=$3

stem=`basename ${file%.fastq.gz}`
report=$report_dir/${stem}_fastqc.zip
if [ -e $report ]; then
    printf "$report already exists.\n"
else
    printf "Generating $report...\n"
    fastqc --outdir $report_dir $file \
	> $log_dir/$stem.out 2> $log_dir/$stem.err
fi

