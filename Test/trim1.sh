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
#   2022-11-17  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 file\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi
trimmed_dir=$1
file=$2

base=`basename $file`
trimmed=$trimmed_dir/$base
if [ -e $trimmed ]; then
    printf "$trimmed already exists.\n"
else
    # Adapter discovered by fastq-scum
    fastq-trim --3p-adapter1 AGATCGGAAGAG --polya-min-length 3 \
	$file $trimmed
fi

