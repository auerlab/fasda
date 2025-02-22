#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for each RNA sample.
##########################################################################

usage()
{
    printf "Usage: $0 replicates\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
replicates=$1

# Document software versions used for publication
uname -a
kallisto version
pwd

gtf=$(Reference/gtf-filename.sh)

# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi

log_dir=Logs/06-kallisto-quant
mkdir -p $log_dir
for r in $(seq 1 $replicates); do
    r2=$(printf "%02d" $r)
    for file in Results/02-trim/cond*-rep$r2.fastq.gz; do
	echo $file
	base=$(basename $file)
	log=$log_dir/${base%.fastq.gz}
	# kallisto 0.46.1 can't handle xz and will simply seg fault rather than
	# issue an error message.  If your trimmed fastq files are in xz format,
	# this will convert to gzip format.
	# Convert xz to gz rather than raw to reduce NFS load from compute nodes
	
	# Kallisto requires an output subdirectory for each sample
	stem=$(basename ${file%.fastq.gz})
	out_dir=Results/06-kallisto-quant/$stem
	mkdir -p $out_dir
    
	if [ ! -e $out_dir/abundance.tsv ]; then
	    set -x
	    kallisto quant \
		--single --fragment-length=190 --sd=10 \
		--threads=$threads \
		--index=Results/05-kallisto-index/transcriptome.index \
		--output-dir=$out_dir $file 2>&1 | tee $log
	    set +x
	else
	    printf "$out_dir/abundance.tsv already exists.\n"
	fi
    done
done
