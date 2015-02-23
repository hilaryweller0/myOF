#!/bin/bash -e

# pre processing
rm -rf 0 0_H
mkdir 0
cp -r system_H system
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0

# set initialisation for H version
setExnerBalancedH

plotPatchData Exner
gv 0/Exner.eps &

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData meshTheta
gv 0/meshTheta.eps &

setUaboveMountain

rm -r system
mv 0 0_H

