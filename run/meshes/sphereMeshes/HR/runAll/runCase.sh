#!/bin/bash -e

base=$1
case=$base

# remove old files and create case
rm -rf $case core
mkdir -p $case
cp -r blankCase/* $case

# create  mesh
sed -i 's/gridN.out/grid'$base'.out/g' $case/system/VoronoiSphereMeshDict
VoronoiSphereMesh -case $case
if [ "$base" -le 5 ] ; then
    plotPatchData -case $case mesh
    gv $case/0/mesh.eps &
fi

# create dual mesh and calculate addressing
polyDualPatch -case $case
meshStructure -case $case '1(originalPatch)'
meshStructure -case $case -region dualMesh '1(originalPatch)'
orthogonality -case $case > $case/orthogLog

# plot primal and dual meshes
if [ "$base" -le 5 ] ; then
    plotPatchData -case $case -time 0 meshUnder
    plotPatchData -case $case -time 0 meshOver -region dualMesh
    cat $case/0/meshOver.eps >> $case/0/meshUnder.eps
    rm $case/0/meshOver.eps
    pstitle $case/0/meshUnder.eps
    makebb $case/0/meshUnder.eps
    gv $case/0/meshUnder.eps &
fi

