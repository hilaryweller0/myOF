#!/bin/bash -e

# running
exnerFoamNonOrthogIMEXEX > log 2>&1 &
tail -f log

# post processing
time=1000
sumFields $time UfDiff $time Uf 0 Uf -scale0 1 -scale1 -1
plotPatchData -time $time thetaUdiff
gv $time/thetaUdiff.eps &


plotPatchData -time $time theta_BF02
gv $time/theta_BF02.eps &
writeuvw -time $time Uf
plotPatchData -time $time w_BF02
gv $time/w_BF02.eps &

sumFields $time thetaDiff $time theta ../../noOrog/dt5_cubicUpFit_1000 theta -scale0 1 -scale1 -1
plotPatchData -time $time thetaDiff
gv $time/thetaDiff.eps &

