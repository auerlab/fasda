#!/bin/sh -e

##########################################################################
#   Description:
#       Compute fold-changes and P-values for N randomly chosen samples
#       from the 48 available for the yeast data.  This demonstrates how
#       variable fold-changes and P-values can be due to biological
#       and technical variation in the coverages.
#       
#   History:
#   Date        Name        Modification
#   2023-05-10  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 sample-count\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
sample_count=$1

# Results/06-kallisto-quant/*/abundance.tsv

