#!/usr/bin/bash

SUFFIX=$1

for FILE in $SUFFIX; do 
  awk 'NR > 1 {if ($3 > max) max = $3} END {print "Max log2FoldChange in", FILENAME ":", max}' "$FILE"
  awk 'NR > 1 {if ($3 < max) max = $3} END {print "Min log2FoldChange in", FILENAME ":", max}' "$FILE"
  echo "_____________________"
done

