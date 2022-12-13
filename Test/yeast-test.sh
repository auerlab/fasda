#!/bin/sh -e

##########################################################################
#   Description:
#       Run all stages of Yeast differential analysis for which
#       output files do not already exist
#       
#   History:
#   Date        Name        Modification
#   2022-11-19  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 max-replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
max_replicates=$1

if [ $max_replicates -lt 3 ]; then
    printf "$0: max-replicates must be at least 3.\n"
    usage
fi

cd ..
./cave-man-install.sh
cd Test

./00-organize.sh
./01-fetch.sh $max_replicates

# FIXME: Let 02-trim.sh handle the count check?
raw_count=$(ls Data/01-fetch/Raw-renamed/*.gz | wc -l)
trimmed_count=$(ls Data/02-trim/*.gz | wc -l)
if [ $trimmed_count -lt $raw_count ]; then
    ./02-trim.sh
else
    printf "02-trim.sh already done.\n"
fi

qc_count=$(ls Data/03-qc/02-trim/*.zip | wc -l)
if [ $qc_count -ne $raw_count ]; then
    ./03-qc.sh
else
    printf "03-qc.sh already done.\n"
fi

if [ ! -e Data/04-reference/all-but-xy.genome.fa.fai ]; then
    ./04-reference.sh
else
    printf "04-reference.sh already done.\n"
fi

if [ ! -e Data/05-kallisto-index/all-but-xy.index ]; then
    ./05-kallisto-index.sh
else
    printf "05-kallisto-index.sh already done.\n"
fi

quant_count=$(ls -d Data/06-kallisto-quant/* | wc -l)
if [ $quant_count -ne $raw_count ]; then
    ./06-kallisto-quant.sh
else
    printf "06-kallisto-quant.sh already done.\n"
fi

./07-fasda-kallisto.sh $max_replicates

if [ ! -e Data/08-hisat-index/all-but-xy.index ]; then
    ./08-hisat-index.sh
else
    printf "08-hisat-index.sh already done.\n"
fi

quant_count=$(ls -d Data/09-hisat-align/* | wc -l)
if [ $quant_count -ne $raw_count ]; then
    ./09-hisat-align.sh
else
    printf "09-hisat-align.sh already done.\n"
fi

./10-fasda-hisat.sh $max_replicates
