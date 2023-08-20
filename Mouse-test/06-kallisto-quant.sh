#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for each RNA sample.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
##########################################################################

# Document software versions used for publication
uname -a
kallisto version
pwd

gtf=$(Reference/gtf-filename.sh)

# If using hdf5, you may need this:
# https://github.com/pachterlab/kallisto/issues/197
# export HDF5_USE_FILE_LOCKING=FALSE

# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi

for zst1 in Results/02-trim/cond*-rep*-R1.fq.zst; do
    zst2=${zst1%-R1.fq.zst}-R2.fq.zst
    echo $zst1 $zst2
    
    # kallisto 0.46.1 can't handle zstd and will simply seg fault rather than
    # issue an error message.
    
    # Kallisto requires an output subdirectory for each sample
    stem=$(basename ${zst1%.fq.zst})
    out_dir=Results/06-kallisto-quant/$stem
    mkdir -p $out_dir

    base1=$(basename $zst1)
    base2=$(basename $zst2)
    
    # kallisto only supports gzip compression as of 0.48.0, so use FIFOs
    # feed it raw fq
    pipe1=/tmp/pipe1-kallisto-$base1
    pipe2=/tmp/pipe2-kallisto-$base2
    rm -f $pipe1 $pipe2
    mkfifo $pipe1 $pipe2
    zstdcat $zst1 > $pipe1 &
    zstdcat $zst2 > $pipe2 &

    set -x
    time kallisto quant \
	--threads=$threads \
	--index=Results/05-kallisto-index/transcriptome.index \
	--output-dir=$out_dir $pipe1 $pipe2
    set +x
    rm -f $pipe1 $pipe2
done
