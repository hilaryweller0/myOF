#!/bin/bash -e

for case in ?/?; do
    nohup nice -n 19 shallowWaterFoam -case $case > $case/log 2>&1 &
    echo $case/log
done

for case in ?/?; do
    cd $case
        ../../../makeRef/makeRef.sh
    cd ../..
done

