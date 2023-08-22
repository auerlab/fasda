#!/bin/sh -e

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
echo Saccharomyces_cerevisiae.R$build.$release.gff3
