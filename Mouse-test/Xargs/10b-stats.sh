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

# Transcripts with modest to high coverage and likely significant
# transcripts='ENSMUST00000000090 ENSMUST00000000109 ENSMUST00000029451'
transcripts=$(sort --random-sort Results/10-fasda-hisat2/fc-1-2-3.txt \
    | awk '$1 != "Feature" && ($2 > 100 || $3 > 100) { print $1 }' | head -n 5)
printf "Randomly selected transcripts:\n$transcripts\n"

for transcript in $transcripts; do
    printf "\nAbundances for $transcript:\n"
    awk -v t=$transcript '$1 == t { split(FILENAME, a, "/"); print a[3], $4 }' \
	Results/10-fasda-hisat2/cond1-rep*
done
pause

ab_dir=Results/10-fasda-hisat2
cd $ab_dir
for transcript in $transcripts; do
    printf "\n=== $transcript ===\n"
    printf "%-20s %10s %10s %10s %10s\n" "Replicate-group" "1" "2" "3" "Mean"
    awk -v t=$transcript \
	'$1 == t { printf("%-20s %10.1f %10.1f %10.1f %10.1f\n", \
		    FILENAME, $2, $3, $4, ($2 + $3 + $4) / 3.0); }' \
	norm-cond1-1-2-3.tsv norm-cond1-4-5-6.tsv
done
pause

for transcript in $transcripts; do
    printf "\n=== $transcript ===\n"
    printf "%-15s %6s %6s %4s %4s %5s %5s %4s\n" \
	"Samples" "MNC1" "MNC2" "SD1" "SD2" "FC" "LFC" "P"
    awk -v t=$transcript \
	'$1 == t { printf("%-15s %6.1f %6.1f %4.1f %4.1f %5.1f %5.1f %4.2f\n",
	    FILENAME, $2, $3, $4, $5, $6, $7, $8); }' \
	*.txt | sort -k 8 -n
done

