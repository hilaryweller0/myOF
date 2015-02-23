#!/bin/bash -e

time=1296000

for case in */[1-9]*; do
    for subCase in $case/save/*diagonalH_CLUSTvort_CLUSTh ; do
        echo post processing case $subCase
        
        if [[ -a $subCase/$time ]]; then
            sumFields -case $subCase $time hError $time h ../../ref/$time h -scale0 1 -scale1 -1
            sumFields -case $subCase $time hTotal $time h constant h0
            ln -sf ../constant/h0 $subCase/$time/h0
            plotPatchData -case $subCase -time $time hError
            plotPatchData -case $subCase -time $time hError200
            #gv $subCase/$time/hError.eps &

            nCells=`grep nFaces $case/constant/polyMesh/boundary | head -1 | \
                    awk '{print $2}' | awk -F';' '{print $1}'`

            ndofs=`grep startFace $case/constant/polyMesh/boundary | tail -1 | \
                    awk '{print $2}' | awk -F';' '{print $1}'`

            dxmax=`grep 'max dx' $case/orthogLog | awk '{print $12}'`
            
            globalSum -case $subCase -time $time hError
            globalSum -case $case/ref -time $time hTotal
            paste $subCase/globalSumhError.dat $case/ref/globalSumhTotal.dat | \
                grep $time | head -1 |\
                awk '{print '$nCells', '$ndofs', $2/$7, $3/$8, $4/$9, '$dxmax'}' \
                    | tee $subCase/hErrors.dat
        fi
    done
done

