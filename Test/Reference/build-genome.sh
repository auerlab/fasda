#!/bin/sh -e

if [ $0 != "Reference/build-genome.sh" ]; then
    cat << EOM

$0 must be run as Reference/build-genome.sh.

EOM
    exit 1
fi

if [ $(uname) = Darwin ]; then
    zcat=gzcat
else
    zcat=zcat
fi

fetch=$(Common/find-fetch.sh)
build=$(Common/genome-build.sh)
release=$(Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# Chromosome files
mkdir -p Data/03-reference
cd Data/03-reference
for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
    file=Saccharomyces_cerevisiae.R$build.dna.chromosome.$chrom.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/saccharomyces_cerevisiae/dna/$file
    fi
done

if [ ! -e $genome ]; then
    printf "Concatenating chromosome FASTAs...\n"
    for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
	printf "$chrom "
	$zcat Saccharomyces_cerevisiae.R$build.dna.chromosome.$chrom.fa.gz >> $genome
    done
    printf "\n"
else
    printf "Using existing $genome...\n"
fi

if [ ! -e $genome.fai ]; then
    printf "Creating index $genome.fai...\n"
    samtools faidx $genome      # Speed up gffread
fi
