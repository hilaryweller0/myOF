# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict constant/polyMesh
blockMesh
add2dMountain
plotPatchData mesh
gv constant/mesh.eps &

# pre processing
rm -rf [0-9]*
mkdir 0
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
plotPatchData meshTheta
gv 0/meshTheta.eps &

mv 0 0_H

# set initialisation for dpdx version
mkdir 0
cp init_0/Exner constant/Exner_init
cp init_0/theta 0/theta
cp init_0/Uf 0
setExnerBalanced
plotPatchData Exner

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

cp 0/theta 0/theta_init
makeHotBubble
plotPatchData meshTheta
gv 0/meshTheta.eps &
mv 0 0_dpdx

# running
exnerFoamNonOrthogIMEXEX > log 2>&1 &
exnerFoamIMEXEX > log 2>&1 &

# post processing
time=1000
writeuvw -time $time Uf
plotPatchData -time $time thetaFilled
gv $time/thetaFilled.eps

sumFields $time thetaDiff $time theta ../../noOrog/1000 theta -scale0 1 -scale1 -1
plotPatchData -time $time thetaDiff
gv $time/thetaDiff.eps &

# make animations
plotPatchData meshTheta -time 0
plotPatchData theta
mkdir -p animations
rm animations/meshTheta*eps
ln -s ../0/meshTheta.eps animations/meshTheta_0.eps
fid=1
for time in ??? ????; do
    if [ -a $time/theta.eps ]; then
        ln -s ../$time/theta.eps animations/meshTheta_$fid.eps
        let fid=$fid+1
    fi
done
