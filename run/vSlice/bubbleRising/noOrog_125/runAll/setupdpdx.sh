#!/bin/bash -e

# pre processing for dpdx initial conditions
rm -rf 0 0_dpdx
mkdir 0
cp -r system_dpdx system
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0

# set initialisation for dpdx version
setExnerBalanced

plotPatchData Exner
gv 0/Exner.eps &

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData meshTheta

rm -r system
mv 0 0_dpdx
gv 0_dpdx/meshTheta.eps &

