#!/bin/bash -e

if [ $# -ne 1 ]; then
    echo usage: runCase caseName
    exit
fi

case=$1

# setup
rm -rf $case
mkdir -p $case
cp -r blankCase$case/* $case
blockMesh -case $case
cp -r $case/constant/polyMesh $case/constant/rMesh

# solve Monge Ampere
MongeAmpere2D -case $case |& tee $case/log

# plot convergence
echo "#time minDet maxDet minEquiDist maxEquiDist detRatio minVol maxVol ratio boostLaplacian" > $case/MAconvergence.dat
grep determinant $case/log | awk '{print $3, $7, $9, $13, $15, $9/$7, $20, $22, $25, $28}'\
    >> $case/MAconvergence.dat
echo "resid" > $case/residsTmp.dat
grep residual $case/log | awk '{print $8}' | awk -F',' '{print $1}'\
     >> $case/residsTmp.dat
paste $case/MAconvergence.dat $case/residsTmp.dat > $case/convAll.dat
rm $case/residsTmp.dat
mv $case/convAll.dat $case/MAconvergence.dat
export case
#gmtPlot plots/MAconvergence.gmt
#gmtPlot plots/MAratios.gmt
gmtPlot plots/boostLaplacian.gmt
gmtPlot plots/MAresids.gmt

# time for other plots
time=`ls $case | sort -n | tail -1`
rm -rf $case/finalTime
mv $case/$time $case/finalTime
rm -r $case/[0-9]*
mv $case/finalTime $case/0
rm -r $case/0/uniform
mv $case/0/rMesh/polyMesh/points $case/constant/polyMesh
rm -r $case/0/rMesh/polyMesh $case/constant/rMesh
mv $case/0/rMesh/* $case/0
rmdir $case/0/rMesh

plotPatchData -case $case -time 0 V
gv $case/0/V.eps &

#plotPatchData -case $case -time 0 source
#gv $case/0/source.eps &

#plotPatchData -case $case -time 0 Phi
#gv $case/0/Phi.eps &

#plotPatchData -case $case -time 0 equidist
#gv $case/0/equidist.eps &

#plotPatchData -case $case -time 0 monitor
#gv $case/0/monitor.eps &

#plotPatchData -case $case -time 0 equidistVol
#gv $case/0/equidistVol.eps &

# mesh analysis
meshAnalysis2D -case $case MongeAmpere2DDict
gmtPlot plots/distDx.gmt
gmtPlot plots/distArea.gmt
gmtPlot plots/distOrthog.gmt
gmtPlot plots/distSkew.gmt

