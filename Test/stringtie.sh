#!/bin/sh

kallisto=Results/06-kallisto-quant/WT-1/abundance.tsv
stringtie=stringtie-transcripts.gtf

if [ ! -e $stringtie ]; then
    stringtie -e \
	-G Results/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
	Results/09-hisat-align/SNF2-1.bam \
	-A stringtie.out \
	-o stringtie.gtf
    awk '$3 == "transcript"' stringtie.gtf > $stringtie
fi

wc -l $kallisto $stringtie
for transcript in $(cat $kallisto | fgrep -v target_id | cut -f 1); do
    echo $transcript
    awk -F ';' -v t="transcript:"$transcript '$2 ~ t { print $3, $4, $5 }' $stringtie
    awk -v t=$transcript '$1 ~ t { print $4 }' $kallisto
    #awk -v t=$transcript '$1 == t { print $8 }' stringtie.out
    echo '==='
done | more
