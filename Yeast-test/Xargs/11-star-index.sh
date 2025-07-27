#!/bin/sh -e

##########################################################################
#   Description:
#       Build star index for reference genome.
##########################################################################

# Document software versions used for publication
uname -a
STAR --version
samtools --version
pwd

reference_dir=Results/04-reference
star_dir=Results/11-star-index

genome=../04-reference/$(Reference/genome-filename.sh)
printf "Using reference $genome...\n"
gtf=../04-reference/$(Reference/gtf-filename.sh)
printf "Using GTF $gtf...\n"

printf "Building STAR index...\n"
mkdir Results/11-star-index
cd Results/11-star-index
STAR \
    --runMode genomeGenerate \
    --genomeSAindexNbases 10 \
    --genomeDir . \
    --genomeFastaFiles $genome \
    --sjdbGTFfile $gtf \
    --sjdbOverhang 49       # Read length - 1
ls
