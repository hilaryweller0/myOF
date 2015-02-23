#!/bin/bash -e

#gmtPlot plots/totalEnergy.gmt
#gmtPlot plots/w.gmt

time=18000
for case in save/*; do
    plotPatchData theta -case $case -time $time
    gv $case/$time/theta.eps &
#    plotPatchData thetaU -case $case -time $time
#    gv $case/$time/thetaU.eps &
    export case=$case
    gmtPlot plots/energy.gmt
    gmtPlot plots/energy2.gmt
    gmtPlot plots/w1.gmt
    gv $case/energy.eps &
    gv $case/energ2.eps &
    gv $case/w.eps &
done

#mkdir -p thetaU
#i=0
#for file in 0/thetaU.eps ???/thetaU.eps ????/thetaU.eps ?????/thetaU.eps; do
#    makebb $file
#    outFile=thetaU/$i.pdf
#    rm -f $outFile
#    #convert -geometry 50%x50% $file $outFile
#    epstopdf --outfile=$outFile $file
#    let i=$i+1
#done

#gv $time/thetaU.eps &

#sumFields $time thetaDiff $time theta 0 theta -scale0 1 -scale1 -1
#plotPatchData -time $time thetaDiff
#gv $time/thetaDiff.eps &

#sumFields $time ExnerDiff $time Exner 0 Exner -scale0 1 -scale1 -1
#plotPatchData -time $time ExnerDiff
#gv $time/ExnerDiff.eps &


