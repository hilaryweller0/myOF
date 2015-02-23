#!/bin/bash -e

schemes=(DubosH_CLUSTPV_CLUSTh
         asymmetricH_CLUSTPV_CLUSTh)

time=432000

gridDirs=(HRbucky/5/save/dt1800
          HRbucky/5_nonOrthog/save/dt1800
          HRbucky/5_Cd/save/dt1800
          diamondCubeC/17x17_eq/save/dt1800
          diamondCubeCd/17x17_eq/save/dt1800
          cubeC/24x24_eq/save/dt1800
          cubeCd/24x24_eq/save/dt1800
)

for scheme in ${schemes[*]}; do
    for case in ${gridDirs[*]}; do
        dir=${case}_${scheme}
        #cd $dir
        #gmtPlot runAll/plotAllDiags
        #cd ../../../..

        # h and U errors
        #sumFields -case $dir $time Udiff $time Uf 0 Uf -scale0 1 -scale1 -1
        sumFields -case $dir $time hDiff $time h 0 h  -scale0 1 -scale1 -1
        plotPatchData -case $dir -time $time hDiff
        gv $dir/$time/hDiff.eps &
    done
done

#for scheme in ${schemes[*]}; do
#    ./comparePlots/errorNorms.sh $scheme
#    ./comparePlots/errorNormsDx.sh $scheme
#done

