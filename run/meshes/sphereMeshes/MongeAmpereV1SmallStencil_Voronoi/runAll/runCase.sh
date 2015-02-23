#!/bin/bash -e

# create Voronoi meshes and their duals from the Monge Ampere meshes

base=$1
res=$2
case=$base/$res
oldDir=MongeAmpereV1SmallStencil

# remove old files and create case
rm -rf $case core
mkdir -p $case
cp -r blankCase/* $case
cp blankCase/constant/requiredResolution ../$oldDir/$case/constant

# create the correct MongeAmpereSphereDict with the correct gamma
gamma=`echo $res | awk '{print (1/$1)**4}'`
sed 's/GAMMA/'$gamma'/g' -i $case/system/MongeAmpereSphereDict

# create a Voronoi mesh and plot
sed -i 's/SOURCECASE/..\/'$oldDir'\/'$base'\/'$res'/g' \
    $case/system/VoronoiSphereMeshDict
VoronoiSphereMesh -case $case |& tee $case/VoronoiMeshLog
#grep Pitteway $case/VoronoiMeshLog > $case/nPitteway.dat

if [ "$base" -le 5 ] ; then
    plotPatchData -case $case -time 0 mesh
    gv $case/0/mesh.eps &
    plotPatchData -case $case -time 0 meshZoom
    gv $case/0/meshZoom.eps &
fi

# create dual mesh and calculate addressing
polyDualPatch -case $case
meshStructure -case $case '1(originalPatch)'
meshStructure -case $case -region dualMesh '1(originalPatch)'
orthogonality -case $case |& tee $case/orthogLog

# mesh analysis
meshAnalysis -case $case MongeAmpereSphereDict
export case
gmtPlot ../plots/distDx.gmt
gmtPlot ../plots/distArea.gmt
gmtPlot ../plots/distOrthog.gmt
gmtPlot ../plots/distCentroidal.gmt
gmtPlot ../plots/distSkew.gmt


