#!/bin/bash -e

# initial conditions
rm -rf 0 0_dpdx*
mv system_dpdx system
mkdir 0
cp init_0/Exner_init constant
cp init_0/theta_init constant
cp init_0/Uf 0
setTheta

# set initialisation for dpdx version
setExnerBalanced

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

mv system system_dpdx
mv 0 0_dpdx


