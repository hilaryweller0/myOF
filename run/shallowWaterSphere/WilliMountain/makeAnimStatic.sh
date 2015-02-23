#!/bin/bash -ve

# post processing for all runs

for time in [0-9]*; do
    ln -sf ../constant/h0 $time/h0
    ln -sf ../../constant/dualMesh/h0 $time/dualMesh/h0
done

vorticityDual

timeStrings=("day_0" "day_0.5" "day_1" "day_1.5" "day_2" "day_2.5"
             "day_3" "day_3.5" "day_4" "day_4.5" "day_5"
             "day_5.5" "day_6" "day_6.5" "day_7" "day_7.5" "day_8"
             "day_8.5" "day_9" "day_9.5" "day_10"
             "day_10.5" "day_11" "day_11.5" "day_12" "day_12.5" "day_13"
             "day_13.5" "day_14" "day_14.5" "day_15")

mkdir -p movie
cp ../../blank/vorticityMovie.tex movie

i=0
times=$(ls -d [0-9]* | sort -n)
for time in $times; do
    sed 's/TIME/'"${timeStrings[$i]}"'/g' ../../../gmtDicts/vorticityTemplate \
        > ../../../gmtDicts/vorticity
    plotPatchData vorticity -time $time -region dualMesh
    convert -density 108 -trim $time/vorticity.eps movie/vorticity$i.jpg
    let i=$i+1
done

cd movie
    pdflatex vorticityMovie.tex
    acroread vorticityMovie.pdf &
cd ..

