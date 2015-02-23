# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp -r system_H system
cp constant/blockMeshDict constant/polyMesh
blockMesh
plotPatchData mesh
gv constant/mesh.eps &

rm -r system

