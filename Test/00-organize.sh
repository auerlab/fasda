#!/bin/sh -e

##########################################################################
#   Description:
#       Create directory structure required for test scripts
#       
#   History:
#   Date        Name        Modification
#   2022-05-12  Jason Bacon Adapt from CNC-EMDiff
##########################################################################

usage()
{
    printf "Usage: $0\n"
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

if [ $# != 0 ]; then
    usage
fi

mkdir -p Data/Raw Data/Raw-renamed Logs
scripts=$(ls 0[1-9]-*) # [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Data/$stage Logs/$stage
done
