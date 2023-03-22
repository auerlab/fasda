#!/bin/sh -e

##########################################################################
#   GFF is used by downstream analysis, such as peak classification
##########################################################################

fetch=$(Common/find-fetch.sh)
release=$(Common/genome-release.sh)
gff=$(Reference/gff-filename.sh)

# macOS zcat looks for .Z extension, while Linux does not have gzcat
zcat='gunzip -c'

##########################################################################
# Ensembl combined GFFs are sorted lexically by chromosome, while BAMs are
# sorted numerically.  Build our own GFF by concatenating individual
# chromosome GFFs in numeric order.  Resorting a GFF is complicated due
# to the hierarchical sort order (all gene components directly under the
# gene, etc).
##########################################################################

cd Results/04-reference
rm -f $gff
site=http://ftp.ensembl.org/pub/release-$release/gff3/saccharomyces_cerevisiae

# Keep header from first GFF.
file=Saccharomyces_cerevisiae.R64-1-1.106.chromosome.I.gff3.gz
if [ ! -e $file ]; then
    printf "Fetching $file...\n"
    $fetch $site/$file
fi
$zcat $file | egrep '^##gff|^#!' | blt deromanize 1 > $gff

# Concatenate the rest without the header
for chrom in II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
    file=Saccharomyces_cerevisiae.R64-1-1.106.chromosome.$chrom.gff3.gz
    if [ ! -e $file ]; then
	printf "Fetching $file...\n"
	$fetch $site/$file
    fi
    $zcat $file | egrep -v '^##[a-z]|^#!' | blt deromanize 1 >> $gff
done

