#!/bin/bash -e

# mesh generation
rm -rf [0-9]*
blockMesh
add2dMountain

# create sponge layer
createSpongeLayer
#plotPatchData sponge
#gv constant/sponge.eps &

