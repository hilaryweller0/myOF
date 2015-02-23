#!/bin/bash -e

gridDirs=(HRbucky/4
          HRbucky/5
          HRbucky/4_nonOrthog
          HRbucky/5_nonOrthog
          cube/12x12_eq
          cube/17x17_eq
          cube/24x24_eq
          cube/48x48_eq
          diamondCube/12x12_eq
          diamondCube/17x17_eq)


for dir in ${gridDirs[*]}; do
    cd $dir
    rm -rf 0
    setWilliSteady
    cd ../..
done




