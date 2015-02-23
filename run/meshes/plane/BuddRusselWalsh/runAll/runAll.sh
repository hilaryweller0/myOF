#!/bin/bash -e

#for res in 2 4 8 16; do
#    ./runAll/runCase.sh X${res}_V1
#done

#for case in Fig*; do
#    ./runAll/runCase.sh $case
#done

#./plots/MAresids.sh

plotPatchData -case Fig4.4_V0 -time 0 Vring
plotPatchData -case Fig4.4_V1 -time 0 Vring
plotPatchData -case Fig4.5_V0 -time 0 Vbell
plotPatchData -case Fig4.5_V1 -time 0 Vbell

for case in Fig*; do
    export case
#    gmtPlot plots/distArea.gmt
#    gmtPlot plots/distDx.gmt
#    gmtPlot plots/distOrthog.gmt
    gmtPlot plots/distSkew.gmt
done

