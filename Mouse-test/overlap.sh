#!/bin/sh -e

file123=Results/10-fasda-hisat2/fc-1-2-3.txt
file124=Results/10-fasda-hisat2/fc-1-2-4.txt
file456=Results/10-fasda-hisat2/fc-4-5-6.txt
fileall=Results/10-fasda-hisat2/fc-all.txt

# Use cat to avoid printing filename.  wc has no option for this.
printf "Total transcripts       = $(cat $file123 | wc -l)\n"
transcripts123=$(awk '$8 < 0.05 { print $1 }' $file123)
printf "P < 0.05 in 1-2-3       = $(echo $transcripts123 | wc -w)\n"
transcripts124=$(awk '$8 < 0.05 { print $1 }' $file124)
printf "P < 0.05 in 1-2-4       = $(echo $transcripts124 | wc -w)\n"
transcripts456=$(awk '$8 < 0.05 { print $1 }' $file456)
printf "P < 0.05 in 4-5-6       = $(echo $transcripts456 | wc -w)\n"
transcriptsall=$(awk '$8 < 0.05 { print $1 }' $fileall)
printf "P < 0.05 in 1-2-3-4-5-6 = $(echo $transcriptsall | wc -w)\n"

awk '$8 < 0.05 { print $1 }' $file456 > sig456
count_both=0
count=0
for transcript in $transcripts123; do
    if fgrep -q $transcript sig456; then
	count_both=$(($count_both + 1))
    fi
    printf "checked: $count  matched: $count_both\r"
    count=$(($count + 1))
done
printf "P < 0.05 in both 1-2-3 and 4-5-6 = $count_both\n"

awk '$8 < 0.05 { print $1 }' $file124 > sig124
count_both=0
count=0
for transcript in $transcripts123; do
    if fgrep -q $transcript sig124; then
	count_both=$(($count_both + 1))
    fi
    printf "checked: $count  matched: $count_both\r"
    count=$(($count + 1))
done
printf "P < 0.05 in both 1-2-3 and 1-2-4 = $count_both\n"
