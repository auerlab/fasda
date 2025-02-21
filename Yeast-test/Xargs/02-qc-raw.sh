#!/bin/sh -e

##########################################################################
#   Description:
#       Trim adapters from raw reads
##########################################################################

usage()
{
    printf "Usage: $0 replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
replicates=$1

# Document software versions used for publication
uname -a
fastqc --version
pwd

# Use fastest compression so gzip can keep up with fastq-trim
export GZIP=-1

hardware_threads=$(./get_hw_threads.sh)
jobs=$(($hardware_threads / 2))
if [ $jobs = 0 ]; then
    jobs=1
fi
printf "Hyperthreads = $hardware_threads  Jobs = $jobs\n"

data_dir=Results/02-qc-raw
log_dir=Logs/02-qc-raw
mkdir -p $data_dir $log_dir
files=""
for r in $(seq 1 $replicates); do
    r2=$(printf "%02d" $r)
    files="$files $(ls Results/01-fetch/Raw-renamed/cond*-rep$r2.fastq.gz)"
done
echo $files | xargs -n 1 -P $jobs ./qc1.sh $data_dir $log_dir
