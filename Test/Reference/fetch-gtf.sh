#!/bin/sh -e

##########################################################################
#   GTF is used by kallisto
##########################################################################

fetch=$(Common/find-fetch.sh)
release=$(Common/genome-release.sh)
gtf=$(Reference/gtf-filename.sh)

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
cd Results/04-reference
if [ ! -e $gtf.gz ]; then
    printf "Fetching $gtf.gz...\n"
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gtf/saccharomyces_cerevisiae/$gtf.gz
fi

if [ ! -e $gtf ]; then
    # Filter for autosomes during decompress
    printf "Uncompressing and filtering $gtf...\n"
    gunzip --stdout $gtf.gz | blt deromanize 1 | awk '$1 ~ "^[0-9]"' > $gtf
fi

