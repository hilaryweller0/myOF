#!/bin/bash -e

time=18000
for case in save/*; do
    if [ -a $case/$time ]; then
        plotPatchData -case $case -time $time w
        gv $case/$time/w.eps &

        sumFields -case $case $time thetaDiff $time theta 0 theta -scale0 1 -scale1 -1
#        sumFields -case $case $time UDiff $time Uf 0 Uf -scale0 1 -scale1 -1
        plotPatchData -case $case -time $time thetaDiff
        gv $case/$time/thetaDiff.eps &
        
        sumFields -case $case $time thetaError $time theta ../../vFine/$time theta -scale0 1 -scale1 -1
        plotPatchData -case $case -time $time thetaError
        gv $case/$time/thetaError.eps &
    fi
done

