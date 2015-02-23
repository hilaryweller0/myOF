#!/bin/bash -e

# initial conditions
rm -rf 0 0_H
mv system_H system
mkdir 0
cp init_0/Exner_init constant
cp init_0/theta_init constant
cp init_0/Uf_init 0/Uf

# set temperature profile
setTheta

# set initialisation for H version
setExnerBalancedH

rm constant/*init
mv system system_H
mv 0 0_H


