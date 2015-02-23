#!/bin/bash -e

time=432000

for num in 3 4 5; do
    for ref in 2 4; do
        case=${num}/${ref}

        if [[ -d $case/$time ]]; then
           echo post processing case $case
            #sumFields -case $case $time Udiff $time Uf 0 Uf -scale0 1 -scale1 -1
            sumFields -case $case $time hDiff $time h 0 h -scale0 1 -scale1 -1
            plotPatchData -case $case -time $time hDiff100
            gv $case/$time/hDiff100.eps &
        fi
    done
done

