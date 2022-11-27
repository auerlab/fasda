#!/bin/sh -e

##########################################################################
#   Description:
#       Run fastqc on raw and trimmed reads
#       
#   History:
#   Date        Name        Modification
#   2022-05-12  Jason Bacon Adapt from CNC-EMDiff
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
