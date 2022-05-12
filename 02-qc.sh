#!/bin/sh -e

report_dir=Data/02-qc/Raw
mkdir -p $report_dir
for file in Data/Yeast/Raw-renamed/*.fastq.gz; do
    stem=`basename ${file%.fastq.gz}`
    report=$report_dir/${stem}_fastqc.zip
    if [ -e $report ]; then
	printf "$report already exists.\n"
    else
	printf "Generating $report...\n"
	fastqc --outdir $report_dir $file
    fi
done
multiqc --outdir $report_dir $report_dir

report_dir=Data/02-qc/Trimmed
mkdir -p $report_dir
for file in Data/01-trim/*.fastq.gz; do
    stem=`basename ${file%.fastq.gz}`
    report=$report_dir/${stem}_fastqc.zip
    if [ -e $report ]; then
	printf "$report already exists.\n"
    else
	printf "Generating $report...\n"
	fastqc --outdir $report_dir $file
    fi
done
multiqc --outdir $report_dir $report_dir
