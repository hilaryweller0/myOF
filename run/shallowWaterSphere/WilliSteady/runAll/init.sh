rm -rf constant/polyMesh constant/dualMesh/polyMesh
mkdir -p constant/dualMesh
ln -s ../../../../../meshes/sphereMeshes/PMAgrid/4_sqrtrho/constant/polyMesh constant/polyMesh
ln -s ../../../../../../meshes/sphereMeshes/PMAgrid/4_sqrtrho/constant/dualMesh/polyMesh constant/dualMesh/polyMesh


rm -rf 0
setWilliSteady
plotPatchData -time 0 hU
gv 0/hU.eps &


