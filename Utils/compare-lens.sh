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

cat << EOM

Note: There will be some mismatches if the kallisto output is not from
the same build/release as the FASDA data.

EOM
pause

while read line; do
    id=$(printf "$line\n" | cut -f 1)
    len=$(printf "$line\n" | cut -f 2)
    kallisto_len=$(awk -v id=$id '$1 == id { print $2 }' \
	Data/Kallisto/chondro-sample1-rep1-time1/abundance.tsv)
    printf "%s %s %s\n" $id $len $kallisto_len
    if [ 0$len != 0$kallisto_len ]; then
	printf "Mismatch.\n"
    fi
done < Data/Hisat2/chondro-sample1-rep1-time1-abundance.tsv

