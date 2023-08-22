#!/bin/sh


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

for condition in WT SNF2; do
    kallisto=Results/06-kallisto-quant/$condition-1/abundance.tsv
    bam=Results/09-hisat2-align/$condition-1.bam
    bam_sorted=${bam%.bam}-sorted.bam
    
    # Make sure hisat2 output is properly sorted
    # Does not seem to be necessary, but is indicated by the stringtie manual
    if [ ! -e $bam_sorted ]; then
	samtools sort -o $bam_sorted $bam
    fi
    
    stringtie_out=$condition-1-stringtie.gtf
    stringtie_transcripts=$condition-1-stringtie-transcripts.gtf
    
    if [ ! -e $stringtie_transcripts ]; then
	stringtie -e -B \
	    -G Results/04-reference/Saccharomyces_cerevisiae.R64-1-1.106.gff3 \
	    $bam_sorted \
	    -o $stringtie_out
	awk '$3 == "transcript"' $stringtie_out > $stringtie_transcripts
    fi
    
    wc -l $kallisto $stringtie_transcripts
done
pause

# Kallisto abundances and stringtie coverages are very different, but
# produce the same fold-changes
reads=$(samtools view Results/09-hisat2-align/SNF2-1.bam | wc -l)
for transcript in $(cat $kallisto | fgrep -v target_id | cut -f 1); do
    echo $transcript
    
    swt=$(awk -F ';' -v t="transcript:"$transcript \
	'$2 ~ t { print $3 }' WT-1-stringtie-transcripts.gtf \
	| awk '{ print $2 }' | tr -d '"')
    
    # Transcript not found
    if [ -z "$swt" ]; then
	continue
    fi
    
    ssnf2=$(awk -F ';' -v t="transcript:"$transcript \
	'$2 ~ t { print $3 }' SNF2-1-stringtie-transcripts.gtf \
	| awk '{ print $2 }' | tr -d '"')

    # Reverse compute raw counts from stringtie FPKM
    start=$(awk -F '\t' -v t="transcript:"$transcript \
	'$9 ~ t { print $4 }' WT-1-stringtie-transcripts.gtf)

    end=$(awk -F '\t' -v t="transcript:"$transcript \
	'$9 ~ t { print $5 }' WT-1-stringtie-transcripts.gtf)
    
    fpkm=$(awk -F ';' -v t="transcript:"$transcript \
	'$2 ~ t { print $4 }' WT-1-stringtie-transcripts.gtf \
	| awk '{ print $2 }' | tr -d '"')
    
    # printf "reads = $reads start = $start end = $end FPKM = $fpkm\n"
    # Just guessing the "/ 2" here.  Single vs paired?
    # wt_count=`printf "$fpkm * ($reads / 1000000) * (($end - $start) / 1000) / 2\n" | bc -l`
    
    # From http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    # regarding prepDE.py3
    #  reads_per_transcript = coverage * transcript_len / read_len
    wt_count=`printf "$swt * ($end - $start) / 50\n" | bc -l`

    start=$(awk -F '\t' -v t="transcript:"$transcript \
	'$9 ~ t { print $4 }' SNF2-1-stringtie-transcripts.gtf)
    end=$(awk -F '\t' -v t="transcript:"$transcript \
	'$9 ~ t { print $5 }' SNF2-1-stringtie-transcripts.gtf)
    fpkm=$(awk -F ';' -v t="transcript:"$transcript \
	'$2 ~ t { print $4 }' SNF2-1-stringtie-transcripts.gtf \
	| awk '{ print $2 }' | tr -d '"')
    # printf "reads = $reads start = $start end = $end FPKM = $fpkm\n"
    
    # Just guessing the "/ 2" here.  Single vs paired?
    # snf2_count=`printf "$fpkm * ($reads / 1000000) * (($end - $start) / 1000) / 2\n" | bc -l`
    
    # From http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    # regarding prepDE.py3
    #  reads_per_transcript = coverage * transcript_len / read_len
    snf2_count=`printf "$ssnf2 * ($end - $start) / 50\n" | bc -l`
    
    # Add a minor error in exchange for avoiding div by 0
    test -z "$swt" && swt=0
    test -z "$ssnf2" && ssnf2=0
    
    sft=`printf "($wt_count + 0.00001) / ($snf2_count + 0.00001)\n" | bc -l 2> /dev/null`
    
    tpm=$(grep $transcript WT-1-stringtie-transcripts.gtf | cut -d ';' -f 5 | awk '{ print $2 }' | tr -d '"')
    
    printf "Tool       %10s %10s %10s %10s\n" WT SNF2 FC "WT TPM"
    
    printf "Stringtie  %10.2f %10.2f %10.2f %10.2f\n" \
	$wt_count $snf2_count $sft $tpm
    
    kwt=$(awk -v t=$transcript '$1 ~ t { print $4 }' \
	Results/06-kallisto-quant/WT-1/abundance.tsv)
    
    ksnf2=$(awk -v t=$transcript '$1 ~ t { print $4 }' \
	Results/06-kallisto-quant/SNF2-1/abundance.tsv)

    wt_tpm=$(awk -v t=$transcript '$1 ~ t { print $5 }' \
	Results/06-kallisto-quant/WT-1/abundance.tsv)

    # Transcript not found
    if [ -z "$kwt" ]; then
	continue
    fi
    
    test -z "$kwt" && kwt=0
    test -z "$ksnf2" && ksnf2=0
    kft=`printf "($kwt + 0.00001) / ($ksnf2 + 0.00001)\n" | bc -l 2> /dev/null`
    printf "Kallisto   %10.2f %10.2f %10.2f %10.2f\n" \
	$kwt $ksnf2 $kft $wt_tpm

    fwt=$(awk -v t=$transcript '$1 ~ t { print $4 }' \
	Results/09-hisat2-align/WT-1-abundance.tsv)
    
    fsnf2=$(awk -v t=$transcript '$1 ~ t { print $4 }' \
	Results/09-hisat2-align/SNF2-1-abundance.tsv)

    printf "FASDA      %10.2f %10.2f\n\n" $fwt $fsnf2
done 2>&1 | more
