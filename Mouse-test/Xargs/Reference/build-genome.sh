#!/bin/sh -e

if [ $0 != "Reference/build-genome.sh" ]; then
    cat << EOM

$0 must be run as Reference/build-genome.sh.

EOM
    exit 1
fi

fetch='curl -O'
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# Chromosome files
mkdir -p Results/04-reference
cd Results/04-reference
for chromosome in $(seq 1 19) X Y; do
    printf "$chrom\n"
    file=Mus_musculus.GRCm$build.dna.chromosome.$chromosome.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna/$file
    fi
done

if [ ! -e $genome ]; then
    printf "Concatenating chromosome FASTAs...\n"
    for chrom in $(seq 1 19) X Y; do
	printf "$chrom "
	# macOS zcat looks for .Z extension, while Linux does not have gzcat
	gunzip -c Mus_musculus.GRCm$build.dna.chromosome.$chrom.fa.gz >> $genome
    done
    printf "\n"
else
    printf "Using existing $genome...\n"
fi

if [ ! -e $genome.fai ]; then
    printf "Creating index $genome.fai...\n"
    samtools faidx $genome      # Speed up gffread
fi
