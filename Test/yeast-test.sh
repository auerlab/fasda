#!/bin/sh -e

cd ..
./cave-man-install.sh
cd Test

./00-organize.sh

raw_count=$(ls Data/Raw-renamed/*.gz | wc -l)
if [ $raw_count -lt 8 ]; then
    printf "How many samples would you like to download? [10] "
    read count
    if [ -z "$count" ]; then
	count=10
    fi
    ./01-fetch.sh $count
fi

trimmed_count=$(ls Data/02-trim/*.gz | wc -l)
if [ $trimmed_count -ne $raw_count ]; then
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

./09-fasda.sh

