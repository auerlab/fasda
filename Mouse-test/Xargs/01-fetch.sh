#!/bin/sh -e

##########################################################################
#   Description:
#       Fetch mouse sample data and create symlinks with descriptive names
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
##########################################################################

usage()
{
    printf "Usage: $0 PRJNA1004253|PRJNA1004652 test|full\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi
prj=$1
mode=$2

case $mode in
test)
    # Just download 1,000,000 spots (fragments)
    # fasterq-dump is newer, but does not support -X clipping
    dump=fastq-dump
    dump_flags='-X 1000000 -v --split-3'
    ;;

full)
    dump=fasterq-dump
    dump_flags='--progress'
    #dump=fastq-dump
    #dump_flags='-v --split-3'
    ;;

*)
    usage
    ;;

esac

# Document software versions used for publication
uname -a
$dump --version || true
pwd

raw=Results/01-fetch/Raw
raw_renamed=Results/01-fetch/Raw-renamed
mkdir -p $raw $raw_renamed

case $prj in
PRJNA1004253)
    conditions="control susceptible"
    ;;

PRJNA1004652)
    conditions="WT CEK"
    ;;

*)
    printf "Invalid PRJ: $prj.\n"
    exit 1
    ;;

esac

cond_num=1
for condition in $conditions; do
    fgrep ",$condition," SraRunTable-$prj.txt | cut -d , -f 1 > $condition.tsv
    printf "$condition:\n"

    biorep=1
    for sample in $(awk '{ print $1 }' $condition.tsv); do
	fq1="${sample}_1.fastq"
	fq2="${sample}_2.fastq"
	printf "$sample = cond$cond_num-rep$biorep...\n"
	if [ ! -e $raw/$fq1.zst ] || [ ! -e $raw/$fq2.zst ]; then
	    if [ ! -e $raw/$fq1 ] || [ ! -e $raw/$fq2 ]; then
		set -x
		# fasterq-dump 3.2.0 will fail if a directory named after the
		# accession does not exist.  mkdir is sufficient, but using
		# prefetch first can speed things up according to the docs.
		# mkdir -p $sample
		prefetch --progress $sample
		$dump $dump_flags --outdir $raw $sample
		set +x
		rm -rf $sample
	    else
		printf "$raw/$fq1 exists.\n"
	    fi
	    
	    # Background so next download can start
	    printf "Compressing...\n"
	    zstd --rm $raw/$fq1 > /dev/null &
	    zstd --rm $raw/$fq2 > /dev/null &
	else
	    printf "$fq1 and $fq2 already exist.\n"
	fi
	(cd $raw_renamed && ln -fs ../Raw/$fq1.zst cond$cond_num-rep$biorep-R1.fq.zst)
	(cd $raw_renamed && ln -fs ../Raw/$fq2.zst cond$cond_num-rep$biorep-R2.fq.zst)
	(cd $raw_renamed && ln -fs cond$cond_num-rep$biorep-R1.fq.zst $condition-rep$biorep-R1.fq.zst)
	(cd $raw_renamed && ln -fs cond$cond_num-rep$biorep-R2.fq.zst $condition-rep$biorep-R2.fq.zst)
	biorep=$(($biorep + 1))
    done
    rm -f $condition.tsv
    cond_num=$(($cond_num + 1))
done
ls -l $raw
ls -l $raw_renamed
