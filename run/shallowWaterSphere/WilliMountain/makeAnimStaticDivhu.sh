#!/bin/bash -ve

# post processing for all runs

for time in [0-9]*; do
    ln -sf ../constant/h0 $time/h0
    ln -sf ../../constant/dualMesh/h0 $time/dualMesh/h0
done

mkdir -p movie

i=0
times=$(ls -d 0 [0-9]??? [0-9]???? | sort -n)
for time in $times; do
    #convert -density 108 -trim $time/divhu.eps movie/divhu$i.jpg
    plotPatchData -time $time -region dualMesh vorticity
    convert -density 108 -trim $time/vorticity.eps movie/vort$i.jpg
    let i=$i+1
done

