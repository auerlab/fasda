#!/bin/sh -e

##########################################################################
#   Description:
#       Build reference genome and transcriptome for all aligners.
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
samtools --version
gffread --version
blt --version
pwd

#############################################################################
# Choose cdna.sh or gtf2fasta.sh, or download the prebuilt kallisto
# file mus_musculus.tar.gz, extract, and copy transcriptome.idx to
# Data/03-kallisto-index/all-but-xy.index.
#
# There are multiple possible transcriptome references that can be used with
# kallisto.
#
# There are many gene IDs referenced in the GTF/GFF that are not in the CDNA.
# FIXME: Document the reason for this.  Are these predicted genes?
#
# There are a few gene IDs in the release 98 CDNA that are not in the GFF.
# If downstream analysis involved looking up genes in the GFF, this could
# result in a few misses.  This caused minor problems for CNC-EMDiff, which
# used CDNA as the reference.
#
# More importantly, the CDNA does not document features such as exons, UTRs,
# etc.  If downstream analysis will examine such features, use the GFF.

# Reference/cdna.sh       # Remove XY from Ensembl cdna transcriptome
Reference/gtf2fasta.sh  # Construct genome minus XY using Ensembl GFF

# Reference/create-chrom-sizes.sh
chrom_sizes="Data/03-reference/chromosome-sizes.tsv"
printf "Generating $chrom_sizes...\n"
blt chrom-lens < Data/03-reference/$(Reference/genome-filename.sh) > $chrom_sizes
cat $chrom_sizes
