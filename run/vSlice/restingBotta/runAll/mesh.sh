#!/bin/bash -e

# mesh generation
mv system_H system
rm -rf 0*
blockMesh
add2dMountain
plotPatchData mesh
gv constant/mesh.eps &
mv system system_H

