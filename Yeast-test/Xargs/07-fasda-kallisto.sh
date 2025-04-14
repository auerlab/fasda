#!/bin/sh -e

##########################################################################
#   Description:
#       Run fasda normalize and fold-change on kallisto abundances
##########################################################################

usage()
{
    printf "Usage: $0 max-replicates\n"
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

output_dir=Results/07-fasda-kallisto
mkdir -p $output_dir
cd $output_dir

# Use fasda built by cave-man-install.sh
PATH=../../../../local/bin:$PATH
export PATH

uname -a
fasda --version
pwd

kallisto_dir=../06-kallisto-quant
samples=$(ls $kallisto_dir | wc -l)
replicates=$(($samples / 2))
log_dir=../../Logs/07-fasda-kallisto
mkdir -p $log_dir

# FIXME: Factor out to fasda-mw.sh?
r0=$(printf '%02d' $replicates)
norm_file=all-norm-$r0.tsv
# Debug rm -f $norm_file
if [ ! -e $norm_file ]; then
    printf "Normalizing condition $condition: $replicates replicates\n"
    files=""
    files="$kallisto_dir/*/abundance.tsv"
    printf "%s\n" $files
    set -x
    time fasda normalize --output $norm_file $files \
	2>&1 | tee $log_dir/normalize-$condition-$r0-MW.out
    set +x
fi
printf "\nCondition $condition normalized counts:\n\n"
head -n 5 $norm_file

printf "\nCondition 1 counts:\n"
cut -f 1-$(($replicates + 1)) $norm_file > cond1-all-norm-$r0.tsv
head -n 5 cond1-all-norm-$r0.tsv

printf "\nCondition 2 counts:\n"
cut -f 1,$(($replicates + 2))-$(($samples + 1)) $norm_file > cond2-all-norm-$r0.tsv
head -n 5 cond2-all-norm-$r0.tsv
pause

de_file=fc-$r0.txt
# Debug rm -f $de_file
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
