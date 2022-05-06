#!/bin/sh -e

##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

./cave-man-install.sh

printf "Hisat2 alignments:\n\n"
data_dir=Data/Hisat2
bams="$data_dir/chondro-sample1-rep1-time1.bam \
    $data_dir/chondro-sample2-rep1-time2.bam \
    $data_dir/chondro-sample3-rep1-time3.bam"

gff=Mus_musculus.GRCm39.106-numeric-sort.gff3
for bam in $bams; do
    printf "\nCalculating abundances for $bam...\n"
    raw_abundance=${bam%.bam}-abundance.tsv
    echo $raw_abundance
    test -e $raw_abundance || time ./abundance "$@" Data/$gff $bam
    
    norm_tsv=${bam%.bam}-norm.tsv
    printf "Normalizing $raw_abundance...\n"
    time ./normalize < $raw_abundance > $norm_tsv
    
    norm_tsvs="$norm_tsvs $norm_tsv"
done
pause

printf "\nComputing fold-change...\n"
time ./fold-change $norm_tsvs | more

printf "\nKallisto alignments:\n\n"

data_dir=Data/Kallisto
norm_tsvs=""
for raw_abundance in $data_dir/*/abundance.tsv; do
    echo $raw_abundance
    norm_tsv=${raw_abundance%.tsv}-norm.tsv
    printf "Normalizing $raw_abundance...\n"
    time ./normalize < $raw_abundance > $norm_tsv
    
    norm_tsvs="$norm_tsvs $norm_tsv"
done

printf "\nComputing fold-change...\n"
time ./fold-change $norm_tsvs | more
