#!/bin/sh -e

cproto -I../local/include *.c > diffanal-protos.h
