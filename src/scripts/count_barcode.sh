#!/bin/bash
gunzip -cf - < $1 | \
    awk -F'\t' '{if ($NF!=""){print $NF}}' | \
    tail -n+2 | \
    sort | \
    uniq -c | \
    awk '{print $2"\t"$1}'
