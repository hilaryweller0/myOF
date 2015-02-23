#!/bin/bash -e

for base in 2 3 4 5 6 7 8; do
    ./runAll/runCase.sh $base
done

