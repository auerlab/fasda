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


header()
{
    printf "\n=============================================================\n"
    printf "$1\n"
    printf "=============================================================\n\n"
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
./cave-man-install.sh || true
cd Test

./00-organize.sh
./01-fetch.sh $max_replicates

# FIXME: Let 02-trim.sh handle the count check?
raw_count=$(ls Data/01-fetch/Raw-renamed/*.gz 2> /dev/null | wc -l)
trimmed_count=$(ls Data/02-trim/*.gz 2> /dev/null | wc -l)
if [ $trimmed_count -lt $raw_count ]; then
    header "Trimming raw reads..."
    ./02-trim.sh
else
    header "02-trim.sh already done."
fi

qc_count=$(ls Data/03-qc/02-trim/*.zip 2> /dev/null | wc -l)
if [ $qc_count -ne $raw_count ]; then
    header "Running FastQC quality checks..."
    ./03-qc.sh
else
    header "03-qc.sh already done."
fi

if [ ! -e Data/04-reference/all-but-xy.genome.fa.fai ]; then
    header "Building reference genome / transcriptome..."
    ./04-reference.sh
else
    header "04-reference.sh already done."
fi

if [ ! -e Data/05-kallisto-index/all-but-xy.index ]; then
    header "Building kallisto index..."
    ./05-kallisto-index.sh
else
    header "05-kallisto-index.sh already done."
fi

quant_count=$(ls -d Data/06-kallisto-quant/* 2> /dev/null | wc -l)
if [ $quant_count -ne $raw_count ]; then
    header "Running kallisto transcriptome alignment / quantification..."
    ./06-kallisto-quant.sh
else
    header "06-kallisto-quant.sh already done."
fi

header "Running FASDA differential analysis on kallisto abundances..."
./07-fasda-kallisto.sh $max_replicates

if [ ! -e Data/08-hisat-index/all-but-xy.index ]; then
    header "Building hisat2 index..."
    ./08-hisat-index.sh
else
    header "08-hisat-index.sh already done."
fi

quant_count=$(ls -d Data/09-hisat-align/* 2> /dev/null | wc -l)
if [ $quant_count -ne $raw_count ]; then
    header "Running hisat2 genome alignment..."
    ./09-hisat-align.sh
else
    header "09-hisat-align.sh already done."
fi

header "Running FASDA differential analysis on hisat2 alignments..."
./10-fasda-hisat.sh $max_replicates
