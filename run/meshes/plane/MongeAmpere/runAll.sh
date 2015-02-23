#!/bin/bash -e

# setup
mkdir -p constant/rMesh
rm -rf [0-9]* constant/rMesh/polyMesh movie/mesh[0-9]*.pdf
blockMesh
cp -r constant/polyMesh constant/rMesh
mkdir 0
cp constant/monitor0 0/monitor
setFields
cp 0/monitor constant
cp constant/Phi 0

# solve Monge Ampere
MongeAmpere2D | tee log

# plot convergence
echo "#time minDet maxDet minEquiDist maxEquiDist detRatio minVol maxVol ratio" > MAconvergence.dat
grep "determinant" log | awk '{print $3, $7, $9, $13, $15, $9/$7, $20, $22, $25}'\
    >> MAconvergence.dat
gmtPlot MAconvergence.gmt
gmtPlot MAratios.gmt

# other plots
time=`ls | sort -n | tail -1`

plotPatchData -time $time detHess
gv $time/detHess.eps&

plotPatchData -time $time Phi
gv $time/Phi.eps&

plotPatchData -time $time -region rMesh monitorUnder
plotPatchData -time $time meshOver
cat $time/meshOver.eps >> $time/monitorUnder.eps
rm $time/meshOver.eps
mv $time/monitorUnder.eps $time/monitorR.eps
makebb $time/monitorR.eps
pstitle $time/monitorR.eps
gv $time/monitorR.eps &

plotPatchData -time 0 monitorUnder
plotPatchData -time $time -region rMesh meshOver
cat $time/meshOver.eps >> 0/monitorUnder.eps
rm $time/meshOver.eps
mv 0/monitorUnder.eps $time/monitor.eps
makebb $time/monitor.eps
pstitle $time/monitor.eps
gv $time/monitor.eps&

plotPatchData -time $time -region rMesh V
gv $time/V.eps &

plotPatchData -time $time equidist
gv $time/equidist.eps &

# make a movie of the meshes
plotPatchData mesh -region rMesh
times=`ls -d [0-9]* | sort -n`
i=0
for dir in $times; do
    echo $i $dir
    epstopdf $dir/mesh.eps
    mv $dir/mesh.pdf movie/mesh$i.pdf
    let i=$i+1
done
#cd movie
#    pdflatex meshMovie.tex
#    acroread meshMovie.pdf &
#cd ..

