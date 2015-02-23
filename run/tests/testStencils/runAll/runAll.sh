# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict constant/polyMesh
blockMesh
plotPatchData mesh
gv 0/mesh.eps &

getStencil 40
plotPatchData stencil
gv 0/stencil.eps &

