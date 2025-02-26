#!/bin/sh -e

#lpjs jobs 1
#lpjs processors-per-job 1
#lpjs threads-per-process processors-per-job
# 
#lpjs pmem-per-processor 181MiB
##############################################################################
# Update PATH on a chimeric cluster (multiple operating systems used for
# compute nodes)
#
# The PATH used by the package manager that installed LPJS (/usr/local for
# FreeBSD ports, usually /usr/pkg or /*/pkg for pkgsrc), is automatically
# prepended to the default PATH.  This is overridden by "#lpjs path", so
# if we use it, we must add all directories ourselves.
#
# Add the default non-priveleged pkgsrc prefix used by auto-pkgsrc-setup.
#
# Caution: Different versions of rsync behave differently with respect
# to creating path components at the destination.  Newer rsync requires
# --mkpath while older ones included with macOS and RHEL do not support
# this flag. Set path to use pkgsrc rsync in ~/Pkgsrc/pkg or /*/pkg.
#lpjs path ~/Pkgsrc/pkg/bin:/opt/pkg/bin:/usr/pkg/bin:/usr/local/bin:/usr/bin:/bin

##########################################################################
#   Description:
#       Fetch Yeast sample data and create symlinks with descriptive names
##########################################################################

##########################################################################
#   Main
##########################################################################

replicates=7

# Document software versions used for publication
uname -a
fastq-dump --version || true
# fasterq-dump --version || true
pwd

raw=Results/01-fetch/Raw
raw_renamed=Results/01-fetch/Raw-renamed
mkdir -p $raw $raw_renamed

# Link raw files to WT-rep or SNF-rep to indicate the biological condition
# Link raw files to condX-repYY for easy and consistent scripting
# I usually make cond1 the control (e.g. wild-type) or first time point
sample_num=1
cond_num=1
for condition in WT SNF2; do
    # Select $replicates replicates
    # Get one technical replicate from each biological replicate
    # Col 2 (Lane) indicates technical rep, use samples where Lane = 1
    # Col 3 is SNF2 mutant or WT
    # Col 4 is biological replicate
    awk -v replicates=$replicates -v condition=$condition \
	'$2 == 1 && $3 == condition && $4 <= replicates' \
	ERP004763_sample_mapping.tsv > $condition.tsv
    printf "$condition:\n"

    for sample in $(awk '{ print $1 }' $condition.tsv); do
	fq="$sample.fastq.zst"
	# Use 2 digits for all replicates in filenames for easier viewing
	if [ ! -e $raw/$fq ]; then
	    # Use rsync if possible on local test platforms.  May not have
	    # sra-tools and pulling from coral saves a lot of bandwidth.
	    printf "Downloading $sample = $condition-$biorep = cond$cond_num-rep$biorep...\n"
	    # fasterq-dump --progress --force --outdir $raw $sample
	    # fasterq-dump now fails unless a directory exists with
	    # the same name as the accession.  Docs say it's faster to
	    # use prefetch first anyway.
	    # https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
	    prefetch --progress $sample
	    fasterq-dump --progress --outdir $raw $sample
	    rm -rf $sample  # Remove .sra files
	    printf "Compressing...\n"
	    # Background so next download can start
	    zstd -f --rm $raw/$sample.fastq
	fi
	biorep=$(awk -v sample=$sample '$1 == sample { print $4 }' $condition.tsv)
	biorep_padded=$(printf "%02d" $biorep)
	sample_num_padded=$(printf "%02d" $sample_num)
	(cd $raw_renamed && ln -fs ../Raw/$fq $condition-$biorep_padded.fastq.zst)
	(cd $raw_renamed && ln -fs ../Raw/$fq sample$sample_num_padded-cond$cond_num-rep$biorep_padded.fastq.zst)
	sample_num=$(($sample_num + 1))
    done
    # rm -f $condition.tsv
    cond_num=$(($cond_num + 1))
done
ls -l $raw
ls -l $raw_renamed
