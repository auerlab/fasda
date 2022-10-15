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

mkdir -p Data Logs
scripts=$(ls 0[1-9]-*) # [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Data/$stage Logs/$stage
done

mkdir -p Data/Raw-renamed
cd Data/Raw-renamed
pwd
for dir in ../../Yeast/ERR*; do
    base=$(basename $dir)
    test -e $dir/$base.fastq.gz || true
    fn=`awk -v id=$base '$1 == id { printf("%s-%s.fastq.gz", $3, $4) }' \
	../../ERP004763_sample_mapping.tsv`
    printf "$dir/$base.fastq.gz -> $fn\n"
    ln -sf $dir/$base.fastq.gz $fn
done

if [ "$(ls WT-* | wc -l)" -lt 8 ] \
    || [ "$(ls SNF2-* | wc -l)" -lt 8 ]; then
    cat << EOM

You must download at least 8 biological replicates of WT and 8 biological
replicates of SNF2 RNA-Seq data from

    https://www.ebi.ac.uk/ena/browser/view/PRJEB5348

and place them in the Yeast directory.  They should have path names such as

    Yeast/ERR458493/ERR458493.fastq.gz

Sample names are mapped in ERP004763_sample_mapping.tsv.  You can select
sample names at the link above and download the lot of them as a .zip file.

EOM
    printf "WT-*:   $(ls WT-* | wc -l)\n"
    printf "SNF2-*: $(ls SNF2-* | wc -l)\n"
    
    printf "\nShowing the first 10 samples of each. "
    pause
    
    awk '$3 == "WT"' ERP004763_sample_mapping.tsv | sort -nk 4 -u | head > list
    awk '$3 == "SNF2"' ERP004763_sample_mapping.tsv | sort -nk 4 -u | head >> list
    sort list
    rm list
    exit 1
fi
