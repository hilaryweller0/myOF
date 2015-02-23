#!/bin/bash -e

# mesh generation
mv system_H system
blockMesh
add2dMountain
plotPatchData mesh
gv constant/mesh.eps &
mv system system_H

