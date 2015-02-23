#!/bin/bash -e

# initial conditions
rm -rf 0 [0-9]*[0-9]

mkdir 0
cp init_0/Exner_init constant
cp init_0/theta_init constant
cp init_0/Uf 0

# set initialisation for H version
isothermalBalance_H

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

rm constant/*init

#thermoVars
#for var in Exner theta BruntFreq; do
#    plotPatchData $var
#    gv 0/$var.eps &
#done

