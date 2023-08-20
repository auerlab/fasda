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
    printf "Usage: $0 data-dir log-dir file\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
trimmed_dir=$1
log_dir=$2
raw1=$3

base1=`basename $raw1`
trimmed1=$trimmed_dir/$base1

raw2=${raw1%-R1.fq.zst}-R2.fq.zst
base2=`basename $raw2`
trimmed2=$trimmed_dir/$base2

if [ -e $trimmed1 ]; then
    printf "$trimmed1 already exists.\n"
else
    base=${base1%-R1.fq.zst}
    # Adapter discovered by fastq-scum
    set -x
    fastq-trim --3p-adapter1 AGATCGGAAGAG --polya-min-length 3 \
	$raw1 $trimmed1 $raw2 $trimmed2 \
	> $log_dir/$base.out 2> $log_dir/$base.err
fi

