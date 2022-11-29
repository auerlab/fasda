#!/bin/sh -e

##########################################################################
#   Description:
#       Run fastqc on raw and trimmed reads
#       
#   History:
#   Date        Name        Modification
#   2022-05-12  Jason Bacon Adapt from CNC-EMDiff
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

uname -a
fastqc --version
pwd

# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi
jobs=$threads
printf "Hyperthreads = $threads  Jobs = $jobs\n"

for dir in Raw-renamed 02-trim; do
    report_dir=Data/03-qc/$dir
    log_dir=Logs/03-qc/$dir
    mkdir -p $report_dir $log_dir
    # Run fastqc jobs in parallel
    ls Data/$dir/*.fastq.gz \
	| xargs -n 1 -P $jobs ./qc1.sh $report_dir $log_dir $file
    if which multiqc; then
	multiqc --outdir $report_dir $report_dir
    fi
done
