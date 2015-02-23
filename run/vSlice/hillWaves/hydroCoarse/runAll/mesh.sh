#!/bin/bash -e

# mesh generation
blockMesh
add2dMountain
plotPatchData mesh
gv constant/mesh.eps &

# create sponge layer
createSpongeLayer
plotPatchData sponge
gv constant/sponge.eps &

