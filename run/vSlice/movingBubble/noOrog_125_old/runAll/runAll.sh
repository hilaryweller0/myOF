# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp -r system_H system
cp constant/blockMeshDict constant/polyMesh
blockMesh
rm -r system

# pre processing for H initial conditions
rm -rf 0 0_H
mkdir 0
cp -r system_H system
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0

# set initialisation for H version
setExnerBalancedH

plotPatchData Exner

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData theta
gv 0/theta.eps &

rm -r system
mv 0 0_H

# pre processing for dpdx initial conditions
rm -rf 0 0_dpdx
mkdir 0
cp -r system_dpdx system
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0

# set initialisation for H version
setExnerBalanced

plotPatchData Exner
gv 0/Exner.eps &

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData meshTheta
gv 0/meshTheta.eps &

rm -r system
mv 0 0_dpdx



# running
exnerFoamNonOrthogIMEXEX > log 2>&1 &

# post processing
time=1000
#sumFields $time thetaDiff $time theta 0 theta_init -scale0 1 -scale1 -1
#plotPatchData -time $time theta
plotPatchData -time $time theta_BF02
gv $time/theta_BF02.eps &
writeuvw -time $time Uf
plotPatchData -time $time w_BF02
gv $time/w_BF02.eps &

sumFields $time UfDiff $time Uf 0 Uf -scale0 1 -scale1 -1
plotPatchData -time $time thetaUdiff
gv $time/thetaUdiff.eps &

plotPatchData -time $time thetaw
gv $time/thetaw.eps &

makecpt -Cjet -T300.001/302.001/0.1 > constant/gmtDicts/white_jet.cpt
plotPatchData -time $time thetaNNS
gv $time/thetaNNS.eps &

plotPatchData -time $time thetaFilled
gv $time/thetaFilled.eps &

gmtPlot constant/gmtPlots/plotw.gmt



