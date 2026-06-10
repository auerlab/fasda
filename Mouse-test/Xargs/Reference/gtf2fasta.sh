#!/bin/sh -e

proper_name="Reference/gtf2fasta.sh"
if [ $0 != $proper_name ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

Reference/fetch-gtf.sh

fetch='curl -O'
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)

transcriptome=$(Reference/transcriptome-filename.sh)
gtf=$(Reference/gtf-filename.sh)
genome=$(Reference/genome-filename.sh)

Reference/build-genome.sh

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Kallisto
cd Results/04-reference
if [ ! -e $transcriptome ]; then
    # gtf_to_fasta is part of tophat, which is obsolete
    # gtf_to_fasta $gtf $genome $transcriptome
    # Maybe?
    # bedtools getfasta -fi $genome -bed $gtf > $transcriptome
    
    # Recommended by Biostar RNA-Seq by Example
    # Warning: couldn't find fasta record for 'MT'!
    # Error: no genomic sequence available (check -g option!).
    # Abort trap (core dumped)
    
    # Generate FASTA containing spliced exons for each transcript from GTF
    printf "Converting $gtf to $transcriptome...\n"
    gffread -w $transcriptome -g $genome $gtf
else
    printf "Using existing $transcriptome...\n"
fi
