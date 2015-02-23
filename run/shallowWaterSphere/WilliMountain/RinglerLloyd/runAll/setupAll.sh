#!/bin/bash -ve

DTs=(1800 900 450 225)
base=RinglerLloyd

i=0
for num in 4 5 6 7; do
    for ref in 1 2 4; do
        case=${num}/${ref}
        rm -rf $case
        mkdir -p $case
        cp -r blank/* $case
        
        sed -i 's/DELTAT/'${DTs[$i]}'/g' $case/system/controlDict

        # link with the mesh
        ln -s ../../../../../../meshes/sphereMeshes/$base/$case/0/polyMesh \
              $case/constant/polyMesh
        ln -s ../../../../../../../meshes/sphereMeshes/$base/$case/0/dualMesh/polyMesh \
              $case/constant/dualMesh/polyMesh
        ln -s ../../../../../meshes/sphereMeshes/$base/$case/orthogLog $case/orthogLog

        # create initial conditions (and plot if on a coarse mesh)
        setWilliMountain -case $case
        setWilliMountain -case $case -region dualMesh
        setWilliSteady -case $case
        ln -s ../constant/h0 $case/0/h0
        ln -s ../../constant/dualMesh/h0 $case/0/dualMesh/h0
        if [ "$num" -le 5 ] ; then
            plotPatchData -case $case -time 0 hU
            gv $case/0/hU.eps &
            plotPatchData -case $case -time 0 h0
            gv $case/0/h0.eps &
        fi
    done
    let i=$i+1
done

