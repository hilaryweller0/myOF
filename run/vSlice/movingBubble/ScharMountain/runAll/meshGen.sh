# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict constant/polyMesh
blockMesh
add2dMountain
checkMesh -constant
plotPatchData mesh
gv constant/mesh.eps &

