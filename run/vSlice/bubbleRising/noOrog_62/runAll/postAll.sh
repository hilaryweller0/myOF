#!/bin/bash -e

# comparisons between time-steps

baseCase=cubicUpCPCFit_H_CN
#baseCase=cubicUpCPCFit_H_RKtheta
dtmin=02
dts=(05 1 2 5 10)
dirRef=dt${dtmin}_${baseCase}

time=1000

it=0
while [ "$it" -ne ${#dts[*]} ]; do
    dt=${dts[$it]}
    dir=save/dt${dt}_${baseCase}
    sumFields -case $dir $time thetaError $time theta ../$dirRef/$time theta \
        -scale0 1 -scale1 -1
    globalSum -case $dir -time $time thetaError
    
    #plotPatchData -case $dir -time $time thetaError
    #gv $dir/$time/thetaError.eps &
    
    let it=$it+1
done

head -1 $dir/globalSumthetaError.dat
it=0
while [ "$it" -ne ${#dts[*]} ]; do
    dt=${dts[$it]}
    dir=save/dt${dt}_${baseCase}

    echo $dt `tail -1 $dir/globalSumthetaError.dat`

    let it=$it+1
done

