#!/bin/sh -e

##########################################################################
#   Description:
#       Run fastqc on raw and trimmed reads
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
replicates=$1

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
    files=""
    for r in $(seq 1 $replicates); do
	r2=$(printf "%02d" $r)
	files="$files $(ls Results/$dir/cond*-rep$r2.fastq.gz)"
    done
    echo $files | xargs -n 1 -P $jobs ./qc1.sh $report_dir $log_dir $file
    if which multiqc && [ ! -e $report_dir/multiqc_data ]; then
	export LC_ALL=en_us-UTF-8
	export LANG=en_us-UTF-8
	multiqc --outdir $report_dir $report_dir
    fi
done
