#!/bin/sh -e


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

# Abundances
kallisto=Results/06-kallisto-quant/cond1-rep01/abundance.tsv
hisat2=Results/09-hisat2-align/cond1-rep01-abundance.tsv

for transcript in $(awk '{ print $1 }' $kallisto); do
    printf '===\n'
    printf "Kallisto: "
    grep $transcript $kallisto
    printf "Hisat2:   "
    grep $transcript $hisat2
done | more

# Normalized
kallisto=Results/07-fasda-kallisto/cond2-all-norm-03.tsv
hisat2=Results/10-fasda-hisat2/cond2-all-norm-03.tsv

for transcript in $(awk '{ print $1 }' $kallisto); do
    printf '===\n'
    printf "Kallisto: "
    grep $transcript $kallisto
    printf "Hisat2:   "
    grep $transcript $hisat2
    k=$(awk -v t=$transcript '$1 == t { print $2 + $3 + $4 }' $kallisto)
    h=$(awk -v t=$transcript '$1 == t { print $2 + $3 + $4 }' $hisat2)
    printf "$h / ($k + .0000001)\n" | bc -l
done | more

# Fold-changes
kallisto=Results/07-fasda-kallisto/fc-03.txt
hisat2=Results/10-fasda-hisat2/fc-03.txt
for transcript in $(awk '{ print $1 }' $kallisto | head -20); do
    echo $transcript
    grep $transcript $kallisto
    grep $transcript $hisat2
done
