#!/bin/bash -e

# calculate mesh diagnostics

#for case in RinglerLloyd/[4-7]/*; do
    #meshAnalysis -case $case VoronoiSphereMeshDict
    
    #export case
    #gmtPlot plots/distDx.gmt
    #gmtPlot plots/distArea.gmt
    #gmtPlot plots/distCentroidal.gmt
    #gmtPlot plots/distOrthog.gmt
    #gmtPlot plots/distSkew.gmt
#done

for base in MongeAmpere*/5 RinglerLloyd/5; do
    for res in 2 4 8 16; do
        case=$base/$res

        #meshAnalysis -case $case MongeAmpereSphereDict
    
        export case
#        gmtPlot plots/distArea.gmt
#       gmtPlot plots/distDx.gmt
#        gmtPlot plots/distHessian.gmt
#        gmtPlot plots/distCentroidal.gmt
        gmtPlot plots/distOrthog.gmt
#        gmtPlot plots/distSkew.gmt
    done
done

