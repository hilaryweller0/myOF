#!/bin/bash -e

# make reference solutions for all cases

#cases="PMA/[1-9]/[1-9] RinglerLloyd/[1-9]/[1-9]"
cases="Walko/[1-9]/[1-9]"

for case in $cases; do
    cd $case
    pwd
    ../../../makeRef/makeRef.sh
    cd ../../..
done

