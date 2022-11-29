#!/bin/sh -e

##########################################################################
#   Description:
#       View technical variance (same biological replicate) in
#       yeast samples.
#       
#   History:
#   Date        Name        Modification
#   2022-11-29  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 \n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

dir=Data/Tech-reps
mkdir -p $dir
samples=$(awk '$3 == "WT" && $4 == 1 { print $1 }' ERP004763_sample_mapping.tsv)
for sample in $samples; do
    if [ ! -e $dir/$sample.fastq.gz ]; then
	if [ ! -e $dir/$sample.fastq ]; then
	    fasterq-dump --progress --force --outdir $dir $sample
	fi
	gzip $dir/$sample.fastq &
    fi
done
wait    # Let all background processes finish

for file in $dir/ERR*[0-9].fastq.gz; do
    trimmed=${file%.fastq.gz}-trimmed.fastq.gz
    if [ ! -e $trimmed ]; then
	fastq-trim --3p-adapter1 AGATCGGAAGAG --polya-min-length 3 \
	    $file $trimmed
    fi
done

if [ ! -e Data/04-reference/all-but-xy.genome.fa.fai ]; then
    ./04-reference.sh
else
    printf "04-reference.sh already done.\n"
fi

if [ ! -e Data/05-kallisto-index/all-but-xy.index ]; then
    ./05-kallisto-index.sh
else
    printf "05-kallisto-index.sh already done.\n"
fi

out_dir=$dir/kallisto-quant
mkdir -p $out_dir
quant_count=$(ls -d $out_dir/* | wc -l)
raw_count=$(ls -d $dir/*-trimmed.fastq.gz | wc -l)
if [ $quant_count -ne $raw_count ]; then
    gtf=$(Reference/gtf-filename.sh)
    threads=2
    for file in $dir/*-trimmed.fastq.gz; do
	stem=${file%-trimmed.fastq.gz}
	sample=$(basename $stem)
	printf "$file\n"
	if [ ! -e $outdir/$sample ]; then
	    set -x
	    kallisto quant \
		--single --fragment-length=190 --sd=10 \
		--genomebam \
		    --gtf=Data/04-reference/$gtf \
		    --chromosomes=Data/04-reference/chromosome-sizes.tsv \
		--threads=$threads \
		--index=Data/05-kallisto-index/all-but-xy.index \
		--output-dir=$out_dir/$sample $file
	    set +x
	else
	    printf "$outdir/$sample already exists.\n"
	fi
    done
else
    printf "kallisto-quant already done.\n"
fi

head -5 $out_dir/ERR45849*/abundance.tsv
