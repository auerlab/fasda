#!/bin/sh -e

uname -a
fastqc --version
pwd

for dir in Raw-renamed 02-trim; do
    report_dir=Data/03-qc/$dir
    mkdir -p $report_dir
    for file in Data/$dir/*.fastq.gz; do
	stem=`basename ${file%.fastq.gz}`
	report=$report_dir/${stem}_fastqc.zip
	if [ -e $report ]; then
	    printf "$report already exists.\n"
	else
	    printf "Generating $report...\n"
	    fastqc --outdir $report_dir $file
	fi
    done
    if which multiqc && which firefox; then
	multiqc --outdir $report_dir $report_dir
	# firefox $report_dir/multiqc_report.html || true
    fi
done
