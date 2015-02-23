#!/bin/bash -e

base=$1
res=$2
case=${base}_anim/$res

# remove old files and create case
rm -rf $case core
mkdir -p $case
cp -r blankCase/* $case

# create initial unit sphere mesh and dual
sed -i 's/gridN.out/grid'$base'.out/g' $case/system/VoronoiSphereMeshDict
VoronoiSphereMesh -case $case
polyDualPatch -case $case
meshStructure -case $case '1(originalPatch)'
meshStructure -case $case -region dualMesh '1(originalPatch)'
#rm -r $case/0/polyMesh/*Centres
sed -i -e 's/empty/patch/g' -e 's/inGroups/\/\/inGroups/g' $case/0/polyMesh/boundary
sed -i -e 's/empty/patch/g' -e 's/inGroups/\/\/inGroups/g' $case/0/dualMesh/polyMesh/boundary

# create the mesh to move (copied from original)
mkdir -p $case/0/rMesh
cp -r $case/0/polyMesh $case/0/rMesh

cp $case/constant/Phi $case/0

# create the correct MongeAmpereSphereDict with the correct gamma
gamma=`echo $res | awk '{print (1/$1)**4}'`
sed 's/GAMMA/'$gamma'/g' -i $case/system/MongeAmpereSphereDict

# solve Monge Ampere
MongeAmpereSphere -case $case |& tee $case/log

# plot the meshes
plotPatchData -case $case -region rMesh mesh

# make pdfs for an animation
mkdir $case/movie
times=$(filename $(ls -d $case/[0-9]*) | sort -n)
i=0
for time in ${times[*]}; do
    echo $case $time $i
    epstopdf $case/$time/mesh.eps
    mv $case/$time/mesh.pdf $case/movie/mesh$i.pdf
    let i=$i+1
done

