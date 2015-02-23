#!/bin/bash -e

# mesh generation
# first generate the ground patch
rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict1 constant/polyMesh/blockMeshDict
blockMesh
cp constant/add2dMountainDict1 constant/add2dMountainDict
add2dMountain
writeMeshObj -patchFaces -constant
mv patchFaces_ground_constant.obj constant/triSurface
rm *obj

# next use the ground patch to generate the snappy mesh
rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict2 constant/polyMesh/blockMeshDict
blockMesh
cp constant/add2dMountainDict2 constant/add2dMountainDict
add2dMountain
snappyHexMesh
plotPatchData meshZoom
gv constant/meshZoom.eps &


# create sponge layer
createSpongeLayer
plotPatchData sponge
gv constant/sponge.eps &
mv system system_H

