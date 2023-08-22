#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
##########################################################################

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

# Document software versions used for publication
uname -a
hisat2 --version
pwd

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# samtools sort dumps temp files in CWD
cd Results/09-hisat2-align

for fastq in ../02-trim/cond*-rep*.fastq.gz; do
    gz1=$(ls $fastq)
    gzb=$(basename $gz1)
    bam=${gzb%.*.*}.bam
    sample=${gzb%.fastq.gz}
    log=../../Logs/09-hisat2-align/$sample.err
    if [ ! -e $bam ]; then
	printf "Running hisat2...\n"
	hisat2 --threads 2 -x ../08-hisat2-index/$genome -U $gz1 \
	    2> $log | samtools sort > $bam
	cat $log
    else
	printf "$bam already exists.\n"
    fi
    
    printf "Indexing $bam...\n"
    if [ ! -e $bam.bai ]; then
	samtools index $bam
    fi
done
