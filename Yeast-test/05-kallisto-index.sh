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
if [ ! -e $transcriptome.fai ]; then
    printf "Building $transcriptome.fai...\n"
    samtools faidx Results/04-reference/$transcriptome
fi

printf "Building kallisto index...\n"
set -x
kallisto index --index=Results/05-kallisto-index/transcriptome.index \
    Results/04-reference/$transcriptome
