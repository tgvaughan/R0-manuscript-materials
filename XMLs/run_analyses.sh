#!/bin/bash

pushd ..

job_prefix=R0paper

if ! [ -d Results ]; then
    mkdir Results
fi

bu=36.5
clock=8e-4

for f in XMLs/*.xml; do
    prefix=`basename $f .xml`

    for i in `seq 1 10`; do 
    beast -seed $i \
          -D clockrate=$clock,burate=$bu \
          -statefile Results/$prefix.$i.state \
          -overwrite \
          $f
    done

done

popd
   
