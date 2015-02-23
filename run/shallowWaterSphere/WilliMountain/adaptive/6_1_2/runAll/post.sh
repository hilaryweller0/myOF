#!/bin/bash -ve

# post processing for all runs

for time in [0-9]*/[0-9]*; do
    ln -sf ../constant/h0 $time/h0
    ln -sf ../../constant/dualMesh/h0 $time/dualMesh/h0
done

timeStrings=("day_0" "day_0.5" "day_1" "day_1.5" "day_2" "day_2.5"
             "day_3" "day_3.5" "day_4" "day_4.5" "day_5"
             "day_5.5" "day_6" "day_6.5" "day_7" "day_7.5" "day_8"
             "day_8.5" "day_9" "day_9.5" "day_10"
             "day_10.5" "day_11" "day_11.5" "day_12" "day_12.5" "day_13"
             "day_13.5" "day_14" "day_14.5" "day_15")

mkdir -p movie

sed 's/TIME/'"${timeStrings[0]}"'/g' ../../gmtDicts/vorticityTemplate \
            > ../../gmtDicts/vorticity
plotPatchData -case 0 vorticity -time 0 -region dualMesh
convert -density 108 -trim 0/0/vorticity.eps movie/vorticity0.jpg

i=1
meshes=$(ls -d [0-9]* | sort -n)
for case in $meshes; do
    time=$(filename $(ls -d $case/[0-9]*) | sort -n | tail -1)
    sed 's/TIME/'"${timeStrings[$i]}"'/g' ../../gmtDicts/vorticityTemplate \
        > ../../gmtDicts/vorticity
    plotPatchData -case $case vorticity -time $time -region dualMesh
    convert -density 108 -trim $case/$time/vorticity.eps movie/vorticity$i.jpg
    let i=$i+1
done

cd movie
    pdflatex vorticityMovie.tex
    acroread vorticityMovie.pdf &
cd ..

