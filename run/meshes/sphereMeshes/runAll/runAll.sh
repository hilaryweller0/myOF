#!/bin/bash

for base in 3 4 5 6 7 ; do
#for base in 3 4 ; do
    for res in 2 4 8 16 ; do
        ../runAll/runCase.sh $base $res
    done
done

for case in RinglerLloyd/5/[1-9]*
do
    mv $case/0/polyMesh/cellCentres  $case/0/polyMesh/cellCentresSave
    mv $case/0/polyMesh/faceCentres  $case/0/polyMesh/faceCentresSave
    meshAnalysis -case $case VoronoiSphereMeshDict
    mv $case/0/polyMesh/cellCentresSave  $case/0/polyMesh/cellCentres
    mv $case/0/polyMesh/faceCentresSave  $case/0/polyMesh/faceCentres
done

# post processing
./plots/MAresids.sh
for case in MongeAmpereV*/5/[1-9]* RinglerLloyd/5/[1-9]*
do
    export case
    gmtPlot plots/distArea.gmt
    gmtPlot plots/distDx.gmt
    gmtPlot plots/distOrthog.gmt
    gmtPlot plots/distSkew.gmt
done

for case in MongeAmpereV*/5/[1-9]*; do
    meshAnalysis -case $case MongeAmpereSphereDict
done

for res in 2 4 8 16; do
    for case in MongeAmpereV*/5/$res RinglerLloyd/5/$res; do
        plotPatchData -case $case -time 0 V$res
        gv $case/0/V$res.eps &
    done
done

for case in MongeAmpere*Voronoi/5/* ; do
    mv $case/0/polyMesh/cellCentres  $case/0/polyMesh/cellCentresSave
    mv $case/0/polyMesh/faceCentres  $case/0/polyMesh/faceCentresSave
    meshAnalysis -case $case MongeAmpereSphereDict
    mv $case/0/polyMesh/cellCentresSave  $case/0/polyMesh/cellCentres
    mv $case/0/polyMesh/faceCentresSave  $case/0/polyMesh/faceCentres
    
    export case
    gmtPlot plots/distArea.gmt
    gmtPlot plots/distDx.gmt
    gmtPlot plots/distOrthog.gmt
    gmtPlot plots/distSkew.gmt
done

