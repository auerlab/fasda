#!/bin/sh -e

##########################################################################
#   Description:
#       Trim adapters from raw reads
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
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

# Document software versions used for publication
uname -a
fastq-trim --version
pwd

# Use fastest compression so gzip can keep up with fastq-trim
export GZIP=-1

hardware_threads=$(./get-hw-threads.sh)
jobs=$(($hardware_threads / 2))
if [ $jobs = 0 ]; then
    jobs=1
fi
printf "Hyperthreads = $hardware_threads  Jobs = $jobs\n"

data_dir=Results/02-trim
log_dir=Logs/02-trim
mkdir -p $data_dir $log_dir
ls Results/01-fetch/Raw-renamed/*.fq.zst \
    | xargs -n 1 -P $jobs ./trim1.sh $data_dir $log_dir

