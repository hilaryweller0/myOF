#!/bin/bash -e

for num in 3 4 5 6 7 8; do
#    for type in '' _Cd _nonOrthog; do
    for type in _Cd; do
        case=${num}${type}
        for subCase in $case/save/*; do
            testPerp -case $subCase > $subCase/testPerpLog
            mkdir $subCase/1
            mv $subCase/0/Ufr $subCase/0/Herror $subCase/1
            sumFields -case $subCase 1 UrError 0 Uf 1 Ufr -scale0 1 -scale1 -1
            plotPatchData -case $subCase -time 1 UrError
        done
    done
done

