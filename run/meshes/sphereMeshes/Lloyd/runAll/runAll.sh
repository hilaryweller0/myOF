#!/bin/bash

for base in 3 4 5 6 7 ; do
    for res in 2 4 8 16 ; do
        ./runAll/runCase.sh $base $res
    done
done

