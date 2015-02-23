#!/bin/bash -e

gridDirs=(HRbucky/3
          HRbucky/3_nonOrthog
          HRbucky/3_Cd
          diamondCubeC/6x6_eq
          diamondCubeCd/6x6_eq
          cubeC/6x6_eq
          cubeCd/6x6_eq
)

for case in ${gridDirs[*]}; do
    orthogonality -case $case
    plotPatchData -case $case orthogonality
    plotPatchData -case $case skewness -region dualMesh
    mv $case/0/orthogonality.eps $case/constant/orthogonality.eps
    cat $case/0/skewness.eps >> $case/constant/orthogonality.eps
    rm $case/0/skewness.eps
    #gv $case/constant/orthogonality.eps &
done

