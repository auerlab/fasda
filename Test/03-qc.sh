#!/bin/sh -e

##########################################################################
#   Description:
#       Run fastqc on raw and trimmed reads
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
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

for dir in 01-fetch/Raw-renamed 02-trim; do
    report_dir=Results/03-qc/$dir
    log_dir=Logs/03-qc/$dir
    mkdir -p $report_dir $log_dir
    # Run fastqc jobs in parallel
    ls Results/$dir/*.fastq.gz \
	| xargs -n 1 -P $jobs ./qc1.sh $report_dir $log_dir $file
    if which multiqc; then
	export LC_ALL=en_us-UTF-8
	export LANG=en_us-UTF-8
	multiqc --outdir $report_dir $report_dir
    fi
done
