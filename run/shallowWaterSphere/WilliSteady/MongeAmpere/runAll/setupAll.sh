#!/bin/bash -ve

DTs=(3600 1800 900 450 225)
base=MongeAmpere

i=0
for num in 3 4 5 6 7; do
    for ref in 2 4; do
        case=${num}/${ref}
        rm -rf $case
        mkdir -p $case
        cp -r blank/* $case
        
        sed 's/DELTAT/'${DTs[$i]}'/g' blank/system/controlDict \
            > $case/system/controlDict

        # link with the mesh
        ln -s ../../../../../../meshes/sphereMeshes/$base/$case/0/polyMesh \
              $case/constant/polyMesh
        ln -s ../../../../../../../meshes/sphereMeshes/$base/$case/0/dualMesh/polyMesh \
              $case/constant/dualMesh/polyMesh
        ln -s ../../../../../meshes/sphereMeshes/$base/$case/orthogLog $case/orthogLog

        # create initial conditions (and plot if on a coarse mesh)
        setWilliSteady -case $case
        if [ "$num" -le 4 ] ; then
            plotPatchData -case $case hU
            gv $case/0/hU.eps &
        fi
    done
    let i=$i+1
done

