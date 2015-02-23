#!/bin/bash -e

# mesh generation
rm -rf [0-9]*
mv system_H system
blockMesh
add2dMountain
plotPatchData meshZoom
gv constant/meshZoom.eps &

# create sponge layer
createSpongeLayer
plotPatchData sponge
gv constant/sponge.eps &
mv system system_H

