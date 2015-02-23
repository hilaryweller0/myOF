#!/bin/bash -e

time=1296000

for case in PMA/[1-9]/[1-9] RinglerLloyd/[1-9]/[1-9] Walko/[1-9]/[1-9]; do
    echo post processing case $case
        
    if [[ -a $case/$time ]]; then
        sumFields -case $case $time hError $time h ref/$time h -scale0 1 -scale1 -1
        sumFields -case $case $time hTotal $time h constant h0
        ln -sf ../constant/h0 $case/$time/h0
        plotPatchData -case $case -time $time hError
        plotPatchData -case $case -time $time hError200
        gv $case/$time/hError.eps &

        nCells=`grep nFaces $case/constant/polyMesh/boundary | head -1 | \
                awk '{print $2}' | awk -F';' '{print $1}'`

        ndofs=`grep startFace $case/constant/polyMesh/boundary | tail -1 | \
                awk '{print $2}' | awk -F';' '{print $1}'`
                
        dxmax=`grep 'max dx' $case/orthogLog | awk '{print $12}'`

        globalSum -case $case -time $time hError
        globalSum -case $case/ref -time $time hTotal
        paste $case/globalSumhError.dat $case/ref/globalSumhTotal.dat | \
            grep $time | head -1 |\
            awk '{print '$nCells', '$ndofs', $2/$7, $3/$8, $4/$9, '$dxmax'}' \
                | tee $case/hErrors.dat
    fi
done

