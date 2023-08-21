#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2021-11-24  Jason Bacon Begin
##########################################################################

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

# Document software versions used for publication
uname -a
hisat2 --version
pwd

genome=$(Reference/genome-filename.sh)

# samtools sort dumps temp files in CWD
cd Results/09-hisat2-align

for zst1 in ../02-trim/cond*-rep*-R1.fq.zst; do
    sample=${zst1%-R1.fq.zst}
    zst2=$sample-R2.fq.zst
    zstb=$(basename $zst1)
    bam=${zstb%-R1.*.*}.bam
    printf "\n===\n\nhisat2: $zst1 $zst2\n"
    log=../../Logs/09-hisat2-align/$(basename $sample).err
    if [ ! -e $bam ]; then
	# Hisat2 doesn't understand zstd and performs seeks on the
	# fasta input, so it must be a regular file, not a pipe
	# Raw files are huge, so use gzip to reduce I/O
	gz1=${zst1%.zst}.gz
	gz2=${zst2%.zst}.gz
	set -x
	zstdcat $zst1 | gzip -1 > $gz1 &
	zstdcat $zst2 | gzip -1 > $gz2
	wait    # Let backgrounded recompress finish
	hisat2 --threads 4 -x ../08-hisat2-index/$genome -1 $gz1 -2 $gz2 \
	    2> $log | samtools sort > $bam
	cat $log
	rm -f $gz1 $gz2
	set +x
    else
	printf "$bam already exists.\n"
    fi
    
    if [ ! -e $bam.bai ]; then
	printf "Indexing $bam...\n"
	samtools index $bam
    fi
done
