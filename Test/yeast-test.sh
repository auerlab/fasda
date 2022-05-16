#!/bin/sh -e

cd ..
./cave-man-install.sh
cd Test

./00-organize.sh

raw_count=$(ls Data/Raw-renamed/*.gz | wc -l)

trimmed_count=$(ls Data/01-trim/*.gz | wc -l)
if [ $trimmed_count -ne $raw_count ]; then
    ./01-trim.sh
else
    printf "01-trim.sh already done.\n"
fi

qc_count=$(ls Data/02-qc/01-trim/*.zip | wc -l)
if [ $qc_count -ne $raw_count ]; then
    ./02-qc.sh
else
    printf "02-qc.sh already done.\n"
fi

if [ ! -e Data/03-reference/all-but-xy.genome.fa.fai ]; then
    ./03-reference.sh
else
    printf "03-reference.sh already done.\n"
fi

if [ ! -e Data/04-kallisto-index/all-but-xy.index ]; then
    ./04-kallisto-index.sh
else
    printf "04-kallisto-index.sh already done.\n"
fi

quant_count=$(ls -d Data/05-kallisto-quant/* | wc -l)
if [ $quant_count -ne $raw_count ]; then
    ./05-kallisto-quant.sh
else
    printf "05-kallisto-quant.sh already done.\n"
fi

dir=Data/05-kallisto-quant
for condition in WT SNF2; do
    ../normalize --output $condition-all-norm.tsv $dir/$condition-*/abundance.tsv
done

../fold-change WT-all-norm.tsv SNF2-all-norm.tsv
