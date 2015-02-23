#!/bin/bash -e

gmtPlot plots/energy.gmt
gmtPlot plots/w.gmt

time=18000
for case in */save/*CFC; do
    plotPatchData theta -case $case -time $time
    gv $case/$time/theta.eps &
    export case=$case
    gmtPlot plots/energy1.gmt
    gmtPlot plots/energy2.gmt
    gmtPlot plots/w1.gmt
done

