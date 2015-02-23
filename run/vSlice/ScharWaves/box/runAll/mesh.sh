#!/bin/bash -e

# mesh generation
blockMesh
plotPatchData meshZoom
gv constant/meshZoom.eps &

# create sponge layer
createSpongeLayer
plotPatchData sponge
gv constant/sponge.eps &

