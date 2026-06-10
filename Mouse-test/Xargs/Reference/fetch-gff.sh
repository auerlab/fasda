#!/bin/sh -e

##########################################################################
#   GFF is used by downstream analysis, such as peak classification
##########################################################################

fetch='curl -O'
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
gff=$(Reference/gff-filename.sh)

##########################################################################
# Ensembl combined GFFs are sorted lexically by chromosome, while BAMs are
# sorted numerically.  Build our own GFF by concatenating individual
# chromosome GFFs in numeric order.  Resorting a GFF is complicated due
# to the hierarchical sort order (all gene components directly under the
# gene, etc).
##########################################################################

cd Results/04-reference
rm -f $gff
site=http://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus

# Keep header from first GFF
gff_chr_1=Mus_musculus.GRCm$build.$release.chromosome.1.gff3.gz
if [ ! -e $gff_chr_1 ]; then
    printf "Fetching $gff_chr_1...\n"
    $fetch $site/$gff_chr_1
fi

# macOS zcat looks for .Z extension, while Linux does not have gzcat
printf "Adding $gff_chr_1 to $gff...\n"
gunzip -c $gff_chr_1 > $gff

# Concatenate the rest without the header
for chrom in $(seq 2 19) X Y; do
    gff_chr_N=Mus_musculus.GRCm$build.$release.chromosome.$chrom.gff3.gz
    if [ ! -e $gff_chr_N ]; then
	printf "Fetching $gff_chr_N...\n"
	$fetch $site/$gff_chr_N
    fi
    printf "Adding $gff_chr_N to $gff...\n"
    gunzip -c $gff_chr_N | egrep -v '^##[a-z]|^#!' >> $gff
done
