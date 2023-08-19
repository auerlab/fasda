#!/bin/sh -e

grep exactly Logs/09-hisat2-align/*.err
grep 'Uniquely.*%' Results/12-star-align/*/Log.final.out

