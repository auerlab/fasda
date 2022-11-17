#!/bin/sh -e

proper_name=Reference/cdna.sh
if [ $0 != "$proper_name" ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

if [ $(uname) = Darwin ]; then
    zcat=gzcat
else
    zcat=zcat
fi

# Need GTF for kallisto quant --genomebam in any case
Reference/fetch-gtf.sh

fetch=$(Common/find-fetch.sh)
build=$(Common/genome-build.sh)
release=$(Common/genome-release.sh)
awk=$(Common/find-awk.sh)
transcriptome=$(Reference/transcriptome-filename.sh)

# Can't guarantee this file will always be available.
# You may need to edit this.
cd Data/04-reference
cdna=Saccharomyces_cerevisiae.R$build.cdna.all.fa.gz
if [ ! -e $cdna ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/fasta/saccharomyces_cerevisiae/cdna/$cdna
else
    printf "$cdna already exists.  Remove and rerun to replace.\n"
fi

set -x
$zcat $cdna | $awk -F : -f ../../Reference/keep-autosomes.awk > $transcriptome
