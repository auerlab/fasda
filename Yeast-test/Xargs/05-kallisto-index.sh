#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
##########################################################################

# Document software versions used for publication
uname -a
kallisto version
samtools --version
pwd

transcriptome=$(Reference/transcriptome-filename.sh)
printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e Results/04-reference/$transcriptome.fai ]; then
    printf "Building $transcriptome.fai...\n"
    samtools faidx Results/04-reference/$transcriptome
else
    printf "$transcriptome.fai already exists.\n"
fi

output_dir=Results/05-kallisto-index
mkdir $output_dir
index=$output_dir/transcriptome.index
if [ ! -e $index ]; then
    printf "Building kallisto index...\n"
    set -x
    kallisto index --index=$index Results/04-reference/$transcriptome
    set +x
else
    printf "$index already exists.\n"
fi
