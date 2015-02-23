#!/bin/bash -e

#scheme=asymH_CLUST50PV_CLUSTh_CLUSTH
#scheme=asymH_CLUST50PV_CLUSTh_midH
#scheme=asymH_CLUST50PV_midhH
#scheme=DubosH_CLUST50PV_midhH

#oldScheme=asymH_CLUST50PV_midhH
#scheme=asymH_CLUST50PV_CLUSTh_midH
time=432000

gridDirs=(HRbucky/4/save/dt3600
          HRbucky/5/save/dt1800
          HRbucky/4_nonOrthog/save/dt3600
          HRbucky/5_nonOrthog/save/dt1800
          cube/12x12_eq/save/dt3600
          cube/17x17_eq/save/dt2400
          cube/24x24_eq/save/dt1800
          cube/48x48_eq/save/dt900
          diamondCube/12x12_eq/save/dt2400
          diamondCube/17x17_eq/save/dt1800)

for case in ${gridDirs[*]}; do
#    mkdir ${case}_${scheme}
#    cd ${case}_${oldScheme}
#    cp -r 0 constant system runAll ../../../../${case}_${scheme}
#    cd ../../../../${case}_${scheme}
#    cp ../../../../runAll/system/fv* system
#    cp ../../../../runAll/system/dualMesh/fv* system/dualMesh
    cd ${case}_${scheme}
    shallowWaterFoam > log 2>&1 &
    cd ../../../..
done

# rerun HRbucky/5 cases with smaller time-step
dtOld=3600
dtNew=1800
for scheme in asymH_CLUST50PV_CLUSTh_CLUSTH \
              asymH_CLUST50PV_CLUSTh_midH \
              asymH_CLUST50PV_midhH DubosH_CLUST50PV_midhH; do
    for grid in HRbucky/5 HRbucky/5_nonOrthog; do
        #mkdir $grid/save/dt${dtNew}_${scheme}
        cd $grid/save/dt${dtOld}_${scheme}
        cp -r constant 0 system runAll ../dt${dtNew}_${scheme}
        cd ../dt${dtNew}_${scheme}
        changeWord $dtOld $dtNew system/controlDict
        shallowWaterFoam > log 2>&1 &
        cd ../../../..
    done
done

# post processing

#for case in ${gridDirs[*]}; do
#    dir=${case}_${scheme}
#    cd $dir
#    gmtPlot runAll/plotAllDiags

#    # h and U errors
#    sumFields $time Udiff $time Uf 0 Uf -scale0 1 -scale1 -1
#    sumFields $time hDiff $time h 0 h  -scale0 1 -scale1 -1
#    plotPatchData -time $time hUdiff

#    cd ../../../..
#done

for scheme in asymH_CLUST50PV_CLUSTh_CLUSTH \
              asymH_CLUST50PV_CLUSTh_midH \
              asymH_CLUST50PV_midhH DubosH_CLUST50PV_midhH
do
    ./comparePlots/errorNorms.sh $scheme
    ./comparePlots/errorNormsDx.sh $scheme
done

