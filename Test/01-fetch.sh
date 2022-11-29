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
#       
#   History:
#   Date        Name        Modification
#   2022-11-17  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 sample_count\n"
    printf "Example: $0 10\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
sample_count=$1

if [ $sample_count -lt 3 ]; then
    printf "Sample count must be at least 3 for Mann-Whitney U-test.\n"
    exit 1
fi

# Document software versions used for publication
uname -a
fasterq-dump --version
pwd

raw=Data/Raw
raw_renamed=Data/Raw-renamed
mkdir -p $raw $raw_renamed

for condition in WT SNF2; do
    # Select $sample_count replicates
    # Get one technical replicate from each biological replicate
    # Col 2 (Lane) indicates technical rep, use samples where Lane = 1
    # Col 3 is SNF2 mutant or WT
    # Col 4 is biological replicate
    awk -v sample_count=$sample_count -v condition=$condition \
	'$2 == 1 && $3 == condition && $4 <= sample_count' \
	ERP004763_sample_mapping.tsv > $condition.tsv
    printf "$condition:\n"

    for sample in $(awk '{ print $1 }' $condition.tsv); do
	fq="$sample.fastq.gz"
	biorep=$(awk -v sample=$sample '$1 == sample { print $4 }' $condition.tsv)
	if [ ! -e $raw/$fq ]; then
	    # Use rsync if possible on local test platforms.  May not have
	    # sra-tools and pulling from coral saves a lot of bandwidth.
	    printf "Downloading $sample = $condition-$biorep...\n"
	    coral=$HOME/Coral/Prog/Src/fasda/Test/Data/Raw/$fq
	    if [ -e $coral ]; then
		rsync -av --partial --progress $coral $raw
	    elif hostname | fgrep -q acadix.biz && ! which fasterq-dump; then
		rsync --partial --progress \
		    coral:Prog/Src/fasda/Test/Data/Raw/$sample.fastq.gz $raw
	    else
		fasterq-dump --progress --force --outdir $raw $sample
		printf "Compressing...\n"
		gzip $raw/$sample.fastq
	    fi
	fi
	(cd $raw_renamed && ln -fs ../Raw/$fq $condition-$biorep.fastq.gz)
    done
    rm -f $condition.tsv
done
ls -l $raw
ls -l $raw_renamed
