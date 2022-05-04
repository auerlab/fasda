#!/bin/sh -e

for file in *.c; do
    cproto -I../local/include $file > ${file%.c}-protos.h
done
