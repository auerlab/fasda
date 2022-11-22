#!/bin/sh -e

for file in *.c; do
    if [ $file != fasda.c ]; then
	cproto -I../local/include $file > ${file%.c}-protos.h
    fi
done
