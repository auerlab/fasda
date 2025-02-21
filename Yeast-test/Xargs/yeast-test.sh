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
    printf "Usage: $0 replicates\n"
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
replicates=$1

if [ $replicates -lt 3 ]; then
    printf "$0: replicates must be at least 3.\n"
    usage
fi

cd ..
./cave-man-install.sh || true
cd Yeast-test

./00-organize.sh
header "Fetching yeast read data..."
time ./01-fetch.sh $replicates

# FIXME: Let 02-trim.sh handle the count check?
raw_count=$(ls Results/01-fetch/Raw-renamed/*.gz 2> /dev/null | wc -l)
trimmed_count=$(ls Results/02-trim/*.gz 2> /dev/null | wc -l)
if [ $trimmed_count -lt $raw_count ]; then
    header "Trimming raw reads..."
    time ./02-trim.sh $replicates
else
    header "02-trim.sh already done."
fi

qc_count=$(ls Results/03-qc/02-trim/*.zip 2> /dev/null | wc -l)
if [ $qc_count -ne $raw_count ]; then
    header "Running FastQC quality checks..."
    time ./03-qc.sh $replicates
else
    header "03-qc.sh already done."
fi

if [ ! -e Results/04-reference/all-but-xy.genome.fa.fai ]; then
    header "Building reference genome / transcriptome..."
    time ./04-reference.sh $replicates
else
    header "04-reference.sh already done."
fi

if [ ! -e Results/05-kallisto-index/all-but-xy.index ]; then
    header "Building kallisto index..."
    time ./05-kallisto-index.sh
else
    header "05-kallisto-index.sh already done."
fi

quant_count=$(ls -d Results/06-kallisto-quant/* 2> /dev/null | wc -l)
if [ $quant_count -ne $raw_count ]; then
    header "Running kallisto transcriptome alignment / quantification..."
    time ./06-kallisto-quant.sh $replicates
else
    header "06-kallisto-quant.sh already done."
fi

header "Running FASDA differential analysis on kallisto abundances..."
time ./07-fasda-kallisto.sh $replicates

if [ ! -e Results/08-hisat2-index/all-but-xy.index ]; then
    header "Building hisat2 index..."
    time ./08-hisat2-index.sh
else
    header "08-hisat2-index.sh already done."
fi

quant_count=$(ls -d Results/09-hisat2-align/* 2> /dev/null | wc -l)
if [ $quant_count -ne $raw_count ]; then
    header "Running hisat2 genome alignment..."
    time ./09-hisat2-align.sh $replicates
else
    header "09-hisat2-align.sh already done."
fi

header "Running FASDA differential analysis on hisat2 alignments..."
time ./10-fasda-hisat2.sh $replicates
