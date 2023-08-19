#!/bin/sh -e

##########################################################################
#   Description:
#       Fetch Yeast sample data and create symlinks with descriptive names
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

# Document software versions used for publication
uname -a
fasterq-dump --version || true
pwd

raw=Results/01-fetch/Raw
raw_renamed=Results/01-fetch/Raw-renamed
mkdir -p $raw $raw_renamed

for condition in WT CEK; do
    fgrep ",$condition," SraRunTable.txt | cut -d , -f 1 > $condition.tsv
    printf "$condition:\n"

    biorep=1
    for sample in $(awk '{ print $1 }' $condition.tsv); do
	fq1="${sample}_1.fastq"
	fq2="${sample}_2.fastq"
	printf "$sample = $condition-$biorep...\n"
	if [ ! -e $raw/$fq1.zst ] || [ ! -e $raw/$fq2.zst ]; then
	    if [ ! -e $raw/$fq1 ] || [ ! -e $raw/$fq2 ]; then
		set -x
		fasterq-dump --outdir $raw --progress $sample
		set +x
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
	(cd $raw_renamed && ln -fs ../Raw/$fq1.zst $condition-$biorep-R1.fq.zst)
	(cd $raw_renamed && ln -fs ../Raw/$fq2.zst $condition-$biorep-R2.fq.zst)
	biorep=$(($biorep + 1))
    done
    rm -f $condition.tsv
done
ls -l $raw
ls -l $raw_renamed
