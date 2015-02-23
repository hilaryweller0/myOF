#!/bin/bash -e

base=$1
res=$2
case=$base/$res

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

# plot convergence
echo "#time minDet maxDet minEquiDist maxEquiDist detRatio minVol maxVol ratio boostLaplacian" > $case/MAconvergence.dat
grep determinant $case/log | awk '{print $3, $7, $9, $13, $15, $9/$7, $20, $22, $25, $28}'\
    >> $case/MAconvergence.dat
echo "resid" > $case/residsTmp.dat
grep residual $case/log | awk '{print $8}' | awk -F',' '{print $1}'\
     >> $case/residsTmp.dat
paste $case/MAconvergence.dat $case/residsTmp.dat > $case/convAll.dat
rm $case/residsTmp.dat
mv $case/convAll.dat $case/MAconvergence.dat
export case
gmtPlot ../plots/MAconvergence.gmt
gmtPlot ../plots/MAratios.gmt
gmtPlot  ../plots/boostLaplacian.gmt
gmtPlot  ../plots/MAresids.gmt

# time for other plots
time=`ls $case | sort -n | tail -1`
mv $case/$time $case/finalTime
mv $case/0/rMesh/polyMesh/boundary $case/0/rMesh/polyMesh/faces \
   $case/0/rMesh/polyMesh/neighbour $case/0/rMesh/polyMesh/owner \
   $case/finalTime/rMesh/polyMesh
rm -r $case/[0-9]*
mv $case/finalTime $case/0
rm -r $case/0/uniform
mv $case/0/rMesh/* $case/0
rmdir $case/0/rMesh

plotPatchData -case $case -time 0 mesh
gv $case/0/mesh.eps &

plotPatchData -case $case -time 0 V${res}
gv $case/0/V${res}.eps &

#plotPatchData -case $case -time 0 source
#gv $case/0/source.eps &

#plotPatchData -case $case -time 0 Phi
#gv $case/0/Phi.eps &

#plotPatchData -case $case -time 0 equidist
#gv $case/0/equidist.eps &

# make the boundaries empty again
sed -i -e 's/patch/empty/g' -e 's/\/\/inGroups/inGroups/g' $case/0/polyMesh/boundary
for file in monitorR Phi volRatio detHess equidist gradPhif; do
    sed -i -e 's/zeroGradient/empty/g' -e 's/calculated/empty/g' $case/0/$file
done

# create dual mesh and calculate addressing
polyDualPatch -case $case
meshStructure -case $case '1(originalPatch)'
meshStructure -case $case -region dualMesh '1(originalPatch)'
orthogonality -case $case > $case/orthogLog

## plot primal and dual meshes
#if [ "$base" -le 5 ] ; then
#    plotPatchData -case $case -time 0 meshUnder
#    plotPatchData -case $case -time 0 meshOver -region dualMesh
#    cat $case/0/meshOver.eps >> $case/0/meshUnder.eps
#    rm $case/0/meshOver.eps
#    pstitle $case/0/meshUnder.eps
#    makebb $case/0/meshUnder.eps
#    gv $case/0/meshUnder.eps &
#fi

#plotPatchData -time 0 mesh
#times=`ls -d [0-9]* | sort -n`
#i=0
#for dir in $times; do
#    echo $i $dir
#    epstopdf $dir/mesh.eps
#    mv $dir/mesh.pdf movie/mesh$i.pdf
#    let i=$i+1
#done
#cd movie
#    pdflatex meshMovie.tex
#    acroread meshMovie.pdf &
#cd ..

# mesh analysis
meshAnalysis -case $case MongeAmpereSphereDict
gmtPlot ../plots/distDx.gmt
gmtPlot ../plots/distArea.gmt
#gmtPlot ../plots/distHessian.gmt
#gmtPlot ../plots/distCentroidal.gmt
gmtPlot ../plots/distOrthog.gmt
gmtPlot ../plots/distSkew.gmt

