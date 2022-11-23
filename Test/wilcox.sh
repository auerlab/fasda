#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-11-23  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 yeast-feature replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi
feature=$1
replicates=$(printf "%02s" $2)

awk -v feature=$feature '$1 ~ feature { print $0 }' \
    Data/09-fasda/WT-all-norm-$replicates.tsv | cut -f 2- > c1.tsv

awk -v feature=$feature '$1 ~ feature { print $0 }' \
    Data/09-fasda/SNF2-all-norm-$replicates.tsv | cut -f 2- > c2.tsv

./wilcox.R

awk -v feature=$feature '$1 ~ feature { print $0 }' \
    Data/09-fasda/WT-SNF2-FC-MW-$replicates.txt

