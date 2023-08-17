#!/bin/sh -e

##########################################################################
#   Description:
#       Run star aligner on each RNA sample.
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

##############################################################################
# Align with star, which can handle splice junctions in RNA reads

# Document software versions used for publication
uname -a
STAR --version
pwd

build=$(Common/genome-build.sh)
release=$(Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# samtools sort dumps temp files in CWD
cd Results/12-star-align
pwd

for fastq in $(ls ../02-trim/*.fastq.gz); do
    base=$(basename $fastq)
    sample=${base%.fastq.gz}
    mkdir -p $sample
    cd $sample
    log=../../../Logs/12-star-align/$sample.err
    bam=Aligned.sortedByCoord.out.bam
    if [ ! -e bam ]; then
	printf "Running star...\n"
	STAR \
	    --runThreadN 1 \
	    --genomeDir ../../11-star-index \
	    --readFilesIn ../$fastq \
	    --readFilesCommand 'gunzip -c' \
	    --outSAMtype BAM SortedByCoordinate
    else
	printf "$bam already exists.\n"
    fi
    
    printf "Indexing $bam...\n"
    if [ ! -e $bam.bai ]; then
	samtools index $bam
    fi
    cd ..
done
