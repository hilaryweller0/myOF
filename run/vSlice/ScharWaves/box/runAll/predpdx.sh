#!/bin/bash -e

# initial conditions
rm -rf 0 0_dpdx
#mv system_dpdx system
mkdir 0
cp init_0/Exner_init constant
cp init_0/theta_init constant
cp init_0/Uf_init 0/Uf

# set temperature profile
setTheta

# set initialisation for H version
setExnerBalanced

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

#time=0
#thermoVars -time $time
#for var in BruntFreq Exner theta T; do
#    plotPatchData -time $time $var
#    gv $time/$var.eps &
#done

#rm constant/*init
#mv system system_dpdx
#mv 0 0_dpdx

