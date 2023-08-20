#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2022-05-12  Jason Bacon Adapt from CNC-EMDiff
##########################################################################

# Document software versions used for publication
uname -a
kallisto version
samtools --version
pwd

# Built from genome and GTF
transcriptome=$(Reference/transcriptome-filename.sh)

# Ensembl cdna
# transcriptome=Mus_musculus.GRCm39.cdna.all.fa

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
