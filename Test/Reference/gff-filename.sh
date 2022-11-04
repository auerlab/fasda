#!/bin/sh -e

build=$(Common/genome-build.sh)
release=$(Common/genome-release.sh)
echo Saccharomyces_cerevisiae.R$build.$release.gff3
