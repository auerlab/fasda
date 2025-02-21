#!/bin/sh -e

if [ $0 != "Reference/build-genome.sh" ]; then
    cat << EOM

$0 must be run as Reference/build-genome.sh.

EOM
    exit 1
fi

# macOS zcat looks for .Z extension, while Linux does not have gzcat
zcat='gunzip -c'

fetch='curl -O'
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# Chromosome files
mkdir -p Results/04-reference
cd Results/04-reference
for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
    file=Saccharomyces_cerevisiae.R$build.dna.chromosome.$chrom.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/saccharomyces_cerevisiae/dna/$file
    fi
done

if [ ! -e $genome ]; then
    printf "Concatenating chromosome FASTAs...\n"
    n=1
    for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
	printf "$chrom "
	# Convert Roman chromosome numbers to Arabic
	$zcat Saccharomyces_cerevisiae.R$build.dna.chromosome.$chrom.fa.gz \
	    | sed -e "s|^>$chrom|>$n|" -e "s|:$chrom|:$n|" >> $genome
	n=$((n + 1))
    done
    printf "\n"
else
    printf "Using existing $genome...\n"
fi

if [ ! -e $genome.fai ]; then
    printf "Creating index $genome.fai...\n"
    samtools faidx $genome      # Speed up gffread
fi
