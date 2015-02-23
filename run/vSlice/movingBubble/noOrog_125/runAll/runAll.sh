# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict constant/polyMesh
blockMesh

# pre processing
rm -rf [0-9]*
mkdir 0
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0

# set initialisation for H version
setExnerBalancedH

plotPatchData Exner
gv 0/Exner.eps &

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData meshTheta
gv 0/meshTheta.eps &

# running
exnerFoamNonOrthogIMEXEX > log 2>&1 &
tail -f log

# post processing
time=1000
sumFields $time UfDiff $time Uf 0 Uf -scale0 1 -scale1 -1
plotPatchData -time $time thetaUdiff
gv $time/thetaUdiff.eps &


plotPatchData -time $time theta_BF02
gv $time/theta_BF02.eps &
writeuvw -time $time Uf
plotPatchData -time $time w_BF02
gv $time/w_BF02.eps &

sumFields $time thetaDiff $time theta ../../noOrog/dt5_cubicUpFit_1000 theta -scale0 1 -scale1 -1
plotPatchData -time $time thetaDiff
gv $time/thetaDiff.eps &

