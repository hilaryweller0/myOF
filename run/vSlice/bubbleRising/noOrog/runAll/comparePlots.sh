#!/bin/bash -e

ref=save/dt01_cubicUpCPCFit_H
time=1000

for case in save/dt02_cubicUpCPCFit_H save/dt05_cubicUpCPCFit_H \
            save/dt1_cubicUpCPCFit_H save/dt2_cubicUpCPCFit_H \
            save/dt5_cubicUpCPCFit_H
do
    sumFields -case $case $time thetaTimeError $time theta \
        ../../$ref/$time theta -scale0 1 -scale1 -1
    dt=`grep deltaT $case/system/controlDict | awk '{print $2}' \
        | awk -F';' '{print $1}'`
    globalSum -case $case -time $time theta
    globalSum -case $case -time $time thetaTimeError
    paste $case/globalSumthetaTimeError.dat $case/globalSumtheta.dat | \
                grep $time | head -1 |\
                awk '{print '$dt', $2/$7, $3/$8, $4/$9}' \
                    | tee $case/thetaTimeErrors.dat
done

echo '#dt l2 linf sum' > comparePlots/thetaTimeErrors.dat
cat save/*/thetaTimeErrors.dat >> comparePlots/thetaTimeErrors.dat
gmtPlot runAll/thetaTimeErrors.gmt

