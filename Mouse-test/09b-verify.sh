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

ls -lh Results/09-hisat2-align/*.bam
pause

for file in Logs/09-hisat2-align/*.err; do
    printf "=== $(basename $file) ===\n"
    fgrep 'aligned concordantly' $file
done | more
