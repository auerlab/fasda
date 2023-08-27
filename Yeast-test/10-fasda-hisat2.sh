#!/bin/sh -e

##########################################################################
#   Description:
#       Run fasda normalize and fold-change on hisat2 alignments
##########################################################################

usage()
{
    printf "Usage: $0 replicates\n"
    exit 1
}


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
replicates=$1

cd Results/10-fasda-hisat2

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH

uname -a
fasda --version
pwd

reference_dir=../04-reference
hisat2_dir=../09-hisat2-align
log_dir=../../Logs/10-fasda-hisat2

##########################################################################
#   Compute abundances
##########################################################################

for condition in 1 2; do
    for r in $(seq 1 $replicates); do
	r2=$(printf "%02d" $r)
	bam=cond$condition-rep$r2.bam
	ab=${bam%.bam}-abundance.tsv
	if [ ! -e $ab ]; then
	    printf "Computing abundances for $condition replicate $r2...\n"
	    set -x
	    time fasda abundance --output-dir . 50 \
		$reference_dir/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
		$hisat2_dir/$bam
	    set +x
	    # tr ' ' '-' < ${bam%.bam}-stringtie-genes.tsv | column -t | sort | more
	fi
	
	column -t $ab | head
	wc $ab
    done
done

# FIXME: Factor out to fasda-mw.sh?
r0=$(printf '%02d' $replicates)
for condition in 1 2; do
    norm_file=cond$condition-all-norm-$r0.tsv
    # Debug
    rm -f $norm_file
    if [ ! -e $norm_file ]; then
	printf "Normalizing condition $condition: $replicates replicates\n"
	files=""
	for r in $(seq 1 $replicates); do
	    r2=$(printf "%02d" $r)
	    files="$files cond$condition-rep$r2-abundance.tsv"
	done
	printf "%s\n" $files
	set -x
	time fasda normalize --output $norm_file $files \
	    2>&1 | tee $log_dir/normalize-$condition-$r0-MW.out
	set +x
    fi
    printf "\nCondition $condition normalized counts:\n\n"
    head $norm_file
done
pause

de_file=fc-$r0.txt
# Debug
rm -f $de_file
if [ ! -e $de_file ]; then
    printf "Computing fold-change for $replicates replicates...\n"
    set -x
    time fasda fold-change \
	--output $de_file \
	cond1-all-norm-$r0.tsv cond2-all-norm-$r0.tsv \
	2>&1 | tee $log_dir/fc-$condition-$r0.out
    set +x
fi

pwd
ls
file=fc-$r0.txt
more $file
printf "\n%-25s %10s %10s\n" "File" "Features" "P < 0.05"
printf "%-25s %10s %10s\n" $file: \
	$(cat $file | wc -l) $(awk '$8 < 0.05' $file | wc -l)
