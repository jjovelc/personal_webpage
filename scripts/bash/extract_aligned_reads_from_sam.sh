!#/usr/bin/bash

INFILE=$1

grep -v "^@" $INFILE | awk '($2 != 4) && ($4 > 0)'
