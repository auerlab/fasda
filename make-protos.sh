#!/bin/sh -e

cproto -I../local/include diffanal.c > diffanal-protos.h
