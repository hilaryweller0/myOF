#!/bin/bash -e

mkdir -p ref
rm -f ref/system ref/constant
ln -s ../system ref/system
ln -s ../constant ref/constant
writeCellData -time 0 h
writeCellData -time 0 -region dualMesh pv
sed 1d 0/h.lonLatz \
    | awk '{pi=4*atan2(1,1); print $1*180/pi+180, $2*180/pi}' \
    > ref/ClonLat
sed 1d 0/dualMesh/pv.lonLatz \
    | awk '{pi=4*atan2(1,1); print $1*180/pi+180, $2*180/pi}' \
    > ref/CdlonLat
rm 0/h.lonLatz 0/dualMesh/pv.lonLatz
nCells=`grep nFaces constant/polyMesh/boundary | head -1 | awk '{print $2}' \
        | awk -F';' '{print $1}'`
nDualCells=`grep nFaces constant/dualMesh/polyMesh/boundary | head -1 \
        | awk '{print $2}' | awk -F';' '{print $1}'`

#times=(0 432000 864000 1296000 1728000 2160000 2592000 3024000 3456000 3888000 4320000)
#stimes=(0000 0120 0240 0360 0480 0600 0720 0840 0960 1080 1200)
#times=(0 432000 864000 1296000)
#stimes=(0000 0120 0240 0360)
times=(1296000)
stimes=(0360)
it=0
while [ $it != ${#times[*]} ]
do
    stime=${stimes[$it]}
    time=${times[$it]}
    echo time step $it running sinterp for time $stime $time
    
    sed 's/xxxx/'$stime'/g' ~/f90/stswm/case5_426/sinterp.nml > sinterp.nml
    sinterp < sinterp.nml
    mkdir -p ref/$time
    cat ../../../makeRef/headFoots/hHeader > ref/$time/hTotal
    echo -e $nCells'\n(' >> ref/$time/hTotal
    awk '{ print $3+0 }' outinterp.dat >> ref/$time/hTotal
    cat ../../../makeRef/headFoots/footer >> ref/$time/hTotal
    
    mkdir 999
    sumFields 999 h ref/$time hTotal constant h0 -scale0 1 -scale1 -1
    mv 999/h ref/$time
    rmdir 999

    sed 's/xxxx/'$stime'/g' ~/f90/stswm/case5_426/sinterpd.nml > sinterp.nml
    sinterp < sinterp.nml
    mkdir -p ref/$time/dualMesh
    cat ../../../makeRef/headFoots/vorticityHeader > ref/$time/dualMesh/vorticity
    echo -e $nDualCells'\n(' >> ref/$time/dualMesh/vorticity
    awk '{ print $3+0 }' outinterp.dat >> ref/$time/dualMesh/vorticity
    cat ../../../makeRef/headFoots/footer >> ref/$time/dualMesh/vorticity
    
    let it=$it+1
done
rm outinterp.dat sinterp.nml
calcPV -case ref
globalSum -case ref h
globalSum -case ref pv -region dualMesh

