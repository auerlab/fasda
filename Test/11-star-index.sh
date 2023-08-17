#!/bin/sh -e

##########################################################################
#   Description:
#       Build star index for reference genome.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2023-08-17  Jason Bacon Begin
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
cd Results/11-star-index
STAR \
    --runMode genomeGenerate \
    --genomeSAindexNbases 10 \
    --genomeDir . \
    --genomeFastaFiles $genome \
    --sjdbGTFfile $gtf \
    --sjdbOverhang 49       # Read length - 1
ls
