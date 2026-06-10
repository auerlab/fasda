#!/bin/sh -e

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
echo $build $release

site=https://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/cdna/
fasta=Mus_musculus.GRCm$build.cdna.all.fa

curl -O $site/$fasta.gz
