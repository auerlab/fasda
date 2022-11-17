#!/bin/sh -e

##########################################################################
#   GFF is used by downstream analysis, such as peak classification
##########################################################################

fetch=$(Common/find-fetch.sh)
release=$(Common/genome-release.sh)
gff=$(Reference/gff-filename.sh)

# GFF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
cd Data/04-reference
rm -f $gff
for chrom in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI; do
    site=http://ftp.ensembl.org/pub/release-$release/gff3/saccharomyces_cerevisiae
    file=Saccharomyces_cerevisiae.R64-1-1.106.chromosome.$chrom.gff3.gz
    if [ ! -e $file ]; then
	printf "Fetching $gff.gz...\n"
	$fetch $site/$file
    fi
    if [ $chrom = I ]; then
	zcat $file | egrep '^##gff|^#!' > $gff
    fi
    zcat $file | egrep -v '^##[a-z]|^#!' | blt deromanize 1 >> $gff
done

