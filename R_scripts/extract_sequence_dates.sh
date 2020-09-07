#!/bin/bash

outbreaks=`ls sequences | cut -d/ -f2 | cut -d. -f1`

for outbreak in $outbreaks; do
    grep '^>' ../sequences/$outbreak.masked | cut -d\| -f3 > ../sequence_dates/${outbreak}_dates.txt
done

