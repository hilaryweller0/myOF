#!/bin/bash -e

ref=noOrog_31/save/dt05_cubicUpCPCFit_H
time=1000
Dx=20e3

for dir in noOrog noOrog_250 noOrog_125 noOrog_62
#for dir in noOrog_250 noOrog_125 noOrog_62
do
    mkdir -p $dir/$ref
    rm -rf $dir/$ref/system
    cp -r $dir/system_H $dir/$ref/system
    sed -i 's/startTime       0/startTime       '$time'/g' \
        $dir/$ref/system/controlDict
    ln -sf ../../../constant $dir/$ref/constant
    cp comparePlots/mapFieldsDict $dir/$ref/system
    mapFields -case $dir/$ref -mapMethod cellVolumeWeight -sourceTime $time \
       -fields '(theta)' -consistent $ref
    #plotPatchData -case $dir/$ref -time $time thetaFilled
    #gv $dir/$ref/$time/thetaFilled.eps &
done

for case in noOrog_250/save/dt4_cubicUpCPCFit_H \
            noOrog_125/save/dt2_cubicUpCPCFit_H \
            noOrog_62/save/dt1_cubicUpCPCFit_H; do
    sumFields -case $case $time thetaError $time theta ../../$ref/$time theta \
        -scale0 1 -scale1 -1
    globalSum -case $case -time $time theta
    globalSum -case $case -time $time thetaError
    nCells=`checkMesh -case $case | grep cells: | awk '{print $2}'`
    dx=`echo $nCells | awk '{print '$Dx'/sqrt($1)}'`
    paste $case/globalSumthetaError.dat $case/globalSumtheta.dat | \
                grep $time | head -1 |\
                awk '{print '$nCells', '$dx', $2/$7, $3/$8, $4/$9}' \
                    | tee $case/thetaErrors.dat
done

cat noOrog_62/save/dt1_cubicUpCPCFit_H/thetaErrors.dat \
    noOrog_125/save/dt2_cubicUpCPCFit_H/thetaErrors.dat \
    noOrog_250/save/dt4_cubicUpCPCFit_H/thetaErrors.dat \
    > comparePlots/thetaErrors2.dat

## plot fine mesh result on coarse mesh
#plotPatchData -time 1000 -case noOrog_31/save/dt1_cubicUpCPCFit_H/ thetaUnder
#plotPatchData  -case noOrog_250 meshOver
#cat noOrog_250/constant/meshOver.eps >> noOrog_31/save/dt1_cubicUpCPCFit_H/1000/thetaUnder.eps
#rm noOrog_250/constant/meshOver.eps
#ps2eps noOrog_31/save/dt1_cubicUpCPCFit_H/1000/thetaUnder.eps
#mv noOrog_31/save/dt1_cubicUpCPCFit_H/1000/thetaUnder.eps.eps noOrog_31/save/dt1_cubicUpCPCFit_H/1000/thetaUnder.eps

