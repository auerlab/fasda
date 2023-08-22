#!/bin/sh -e

##########################################################################
#   This script should no longer be needed, since kallisto -gtf appears
#   to work with GFF3.  Use the equivalent gff3 script instead.
##########################################################################

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
echo Saccharomyces_cerevisiae.R$build.$release.gtf
