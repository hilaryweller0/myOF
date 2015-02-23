#!/bin/bash -ve

# monitor funciton minimum value
minVal=1e-5

## create monitor functions for each time-step
#mkdir -p resolutionField
##ncFileIn=/home/hilary/work/pptData/X134.225.102.228.28.5.26.21.nc
#ncFileIn=/home/hilary/work/pptData/X134.225.102.228.28.11.7.47.nc
#t=0
#while [ $t -le 11 ]; do
#    mkdir -p resolutionField/$t
#    ncks -v prate --dimension time,$t,$t $ncFileIn \
#      | grep 'prate\[' | awk '{print $4}' | awk -F'=' '{print $2 + '$minVal'}' \
#      > resolutionField/$t/monitor.dat
#    let t=$t+1
#done

## create lat-lon mesh for the monitor function files
#sphPolarMesh -case resolutionField
#gedit -s resolutionField/constant/polyMesh/boundary

## create the monitor funcitons
#nCells=17666
#t=0
#while [ $t -le 11 ]; do
#    cp resolutionField/constant/monitorHeader resolutionField/$t/monitor
#    echo -e $nCells\\n'('\\n$minVal >> resolutionField/$t/monitor
#    cat resolutionField/$t/monitor.dat >> resolutionField/$t/monitor
#    echo $minVal >> resolutionField/$t/monitor
#    cat resolutionField/constant/monitorFooter >> resolutionField/$t/monitor
#    let t=$t+1
#done


# set up rMeshes for each of the monitor functions
case=.
rm -rf $case/[0-9]*
mkdir $case/0
cp $case/constant/Phi $case/0
cp $case/constant/requiredResolution ../../HR/6/constant/
VoronoiSphereMesh -case $case
polyDualPatch -case $case
meshStructure -case $case '1(originalPatch)'
meshStructure -case $case -region dualMesh '1(originalPatch)'
sed -i -e 's/empty/patch/g' -e 's/inGroups/\/\/inGroups/g' $case/0/polyMesh/boundary
sed -i -e 's/empty/patch/g' -e 's/inGroups/\/\/inGroups/g' $case/0/dualMesh/polyMesh/boundary
# create the mesh to move (copied from original)
mkdir -p $case/0/rMesh
cp -r $case/0/polyMesh $case/0/rMesh
cp $case/constant/Phi $case/0

# create rMeshes for each of the monitor functions
t=0
while [ $t -le 11 ]; do
    rm -rf [1-9]*
    # create monitor field in time directory
    mapFields -fields '(monitor)' -mapMethod cellVolumeWeight -sourceTime $t \
               resolutionField -case $case
    MongeAmpereSphere -case $case

    # create plots save latest time and delete the others
    time=`ls -d [0-9]* | sort -n | tail -1`
    plotPatchData -time $time -region rMesh pptMonitor -case $case
    pscoast -R -J -Wthick,navy -Dc -O >> $case/$time/pptMonitor.eps
    makebb $case/$time/pptMonitor.eps

    mv $case/$time $case/time$t
    gv $case/time$t/pptMonitor.eps &
    rm -rf $case/[1-9]*

    let t=$t+1
done

eps2gif pptMonitor.gif time?/pptMonitor.eps time??/pptMonitor.eps


