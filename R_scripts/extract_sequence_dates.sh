#!/bin/bash

mkdir -p sequence_dates

for seqfile in sequences/*.masked; do
    outbreak=`basename $seqfile .masked`
    grep '^>' sequences/$outbreak.masked | cut -d\| -f3 > sequence_dates/${outbreak}_dates.txt
done

