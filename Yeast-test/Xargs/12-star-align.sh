#!/bin/sh -e

##########################################################################
#   Description:
#       Run star aligner on each RNA sample.
##########################################################################

##############################################################################
# Align with star, which can handle splice junctions in RNA reads

# Document software versions used for publication
uname -a
STAR --version
pwd
which STAR

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# samtools sort dumps temp files in CWD
mkdir Results/12-star-align
cd Results/12-star-align
pwd

for fastq in $(ls ../02-trim/cond*-rep*.fastq.gz); do
    base=$(basename $fastq)
    sample=${base%.fastq.gz}
    mkdir -p $sample
    cd $sample
    log=../../../Logs/12-star-align/$sample.err
    bam=Aligned.sortedByCoord.out.bam
    if [ ! -e bam ]; then
	printf "===\nSTAR $fastq...\n\n"
	STAR \
	    --runThreadN 2 \
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
