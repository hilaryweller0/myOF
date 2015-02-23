#!/bin/bash -e

time=15000
for case in save/*CPC; do
    if [ -a $case/$time ]; then
        plotPatchData -case $case -time $time w
        gv $case/$time/w.eps &

        sumFields -case $case $time thetaDiff $time theta 0 theta -scale0 1 -scale1 -1
        //sumFields -case $case $time UDiff $time Uf 0 Uf -scale0 1 -scale1 -1
        plotPatchData -case $case -time $time thetaDiff
        gv $case/$time/thetaDiff.eps &
    fi
done

