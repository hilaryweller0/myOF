#!/bin/bash -e

base=$1
res=$2
case=$base/$res

# remove old files and create case
rm -rf $case core
mkdir -p $case
cp -r blankCase/* $case

# create initial unit sphere mesh
sed 's/gridN.out/grid'$base'.out/g' $case/system/VoronoiSphereMeshDict1 \
    > $case/system/VoronoiSphereMeshDict
VoronoiSphereMesh -case $case

# create the required resolution
# create the correct MongeAmpereSphereDict with the correct min and max dx
mindx=100000
maxdx=`echo | awk '{print '$res'*'$mindx' }'`
sed -e 's/MAXDX/'$maxdx'/g' -e 's/MINDX/'$mindx'/g' -i $case/system/setResolutionDict
setResolution -case $case

# generate Lloyd algorithm meshes
sed -e 's/BASE/'$base'/g' -e 's/RES/'$res'/g' $case/system/VoronoiSphereMeshDict2 \
    > $case/system/VoronoiSphereMeshDict
VoronoiSphereMesh -case $case

# plot mesh generation process
#times=$(filename $(ls -d $case/*0) | sort -n)
times=$(filename $(ls -d $case/[0-9]*) | sort -n)
i=0
for time in $times; do
    echo $case $i $time
    plotPatchData -case $case -time $time mesh
    epstopdf $case/$time/mesh.eps
    mv $case/$time/mesh.pdf $case/movie/mesh$i.pdf
    let i=$i+1
done
cd $case/movie
    pdflatex meshMovie.tex
    acroread meshMovie.pdf &
cd ../../..


