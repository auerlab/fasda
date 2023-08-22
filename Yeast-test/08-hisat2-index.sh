#!/bin/sh -e

##########################################################################
#   Description:
#       Build hisat2 index for reference genome.
##########################################################################

# Document software versions used for publication
uname -a
hisat2 --version
samtools --version
pwd

# Run hisat2-build on a copy in 08-hisat2-index so it will put the .ht2
# files there
reference_dir=Results/04-reference
hisat2_dir=Results/08-hisat2-index

genome=$(Reference/genome-filename.sh)
ln -f $reference_dir/$genome $hisat2_dir
ln -f $reference_dir/$genome.fai $hisat2_dir
genome=$hisat2_dir/$genome
printf "Using reference $genome...\n"

if [ ! -e $genome.8.ht2 ]; then
    printf "Building $genome.*.ht2...\n"
    hisat2-build $genome $genome
fi
ls $hisat2_dir
