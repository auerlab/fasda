#!/bin/sh -e

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)

echo Mus_musculus.GRCm$build.$release.chr.gtf
