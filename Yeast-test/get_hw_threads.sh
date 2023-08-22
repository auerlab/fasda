#!/bin/sh -e

# FIXME: parallel --number-of-cores reports 2 on unixdev1, should be 16
# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
if [ $(uname) = Linux ]; then
    threads=$(getconf _NPROCESSORS_ONLN)
else
    threads=$(getconf NPROCESSORS_ONLN)
fi
printf "$threads\n"
