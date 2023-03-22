#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       
#   Arguments:
#       
#   Returns:
#
#   Examples:
#
#   Files:
#
#   Environment:
#
#   See also:
#       
#   History:
#   Date        Name        Modification
#   2022-11-22  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 feature-name replicates\n"
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

printf "Replicates = $replicates\n"
awk -v feature=$feature '$1 == feature { print $0 }' \
    Results/09-fasda/WT-all-norm-$replicates.tsv \
    Results/09-fasda/SNF2-all-norm-$replicates.tsv

