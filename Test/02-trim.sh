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

# Document software versions used for publication
uname -a
fastq-trim --version
pwd

# Use fastest compression so gzip can keep up with fastq-trim
export GZIP=-1

# FIXME: parallel --number-of-cores reports 2 on unixdev1, should be 16
# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi
jobs=$(($threads / 2))
if [ $jobs = 0 ]; then
    jobs=1
fi
printf "Hyperthreads = $threads  Jobs = $jobs\n"

data_dir=Data/02-trim
log_dir=Logs/02-trim
mkdir -p $data_dir $log_dir
ls Data/01-fetch/Raw-renamed/*.fastq.gz \
    | xargs -n 1 -P $jobs ./trim1.sh $data_dir $log_dir

