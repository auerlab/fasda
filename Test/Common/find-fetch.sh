#!/bin/sh -e

for fetch in fetch curl wget; do
    if which $fetch > /dev/null 2> /dev/null; then
	break
    fi
done
case $fetch in
curl)
    fetch="curl -O"
    ;;
fetch|wget)
    ;;
*)
    printf "No fetch, curl, or wget found.\n"
    exit 1
    ;;
esac
printf "$fetch\n"
