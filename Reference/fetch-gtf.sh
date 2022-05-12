#!/bin/sh -e

##########################################################################
#   This script should no longer be needed, since kallisto -gtf appears
#   to work with GFF3.  Use the equivalent gff3 script instead.
##########################################################################

fetch=$(Common/find-fetch.sh)
release=$(Common/genome-release.sh)
gtf=$(Reference/gtf-filename.sh)

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
cd Data/03-reference
if [ ! -e $gtf.gz ]; then
    printf "Fetching $gtf.gz...\n"
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gtf/saccharomyces_cerevisiae/$gtf.gz
fi

if [ ! -e $gtf ]; then
    # Filter for autosomes during decompress
    printf "Uncompressing and filtering $gtf...\n"
    gunzip --stdout $gtf.gz | awk '$1 ~ "^[0-9]"' > $gtf
fi
