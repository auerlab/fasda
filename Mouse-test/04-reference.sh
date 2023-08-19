#!/bin/sh -e

##########################################################################
#   Description:
#       Build reference genome and transcriptome for all aligners.
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2019-10-??  Jason Bacon Begin
##########################################################################

# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=4g
#SBATCH --output=Logs/04-reference/slurm-%A.out
#SBATCH --error=Logs/04-reference/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
samtools --version
gffread --version
blt --version
pwd

#############################################################################
# Choose cdna.sh or gtf2fasta.sh, or download the prebuilt kallisto
# file mus_musculus.tar.gz, extract, and copy transcriptome.idx to
# Results/03-kallisto-index/all-but-xy.index.
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
#
# Remove XY from Ensembl cdna transcriptome
# Reference/cdna.sh
#
# Construct transciptome minus XY using GTF and genome
Reference/gtf2fasta.sh

# GFF needed for fasda abundance
Reference/fetch-gff.sh

# Reference/create-chrom-sizes.sh
chrom_sizes="Results/04-reference/chromosome-sizes.tsv"
printf "Generating $chrom_sizes...\n"
blt chrom-lens < Results/04-reference/$(Reference/genome-filename.sh) > $chrom_sizes
cat $chrom_sizes
