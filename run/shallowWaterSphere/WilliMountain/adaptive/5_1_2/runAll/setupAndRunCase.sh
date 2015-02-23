#!/bin/bash -ve

# run the flow over a mountain test case on a variety of differnet meshes, switching
# between meshes and remapping

baseMesh=../../../../../../meshes/sphereMeshes/RinglerLloydThick
meshes=(5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2
        5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2
        5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2 5/1 5/2)
BIGDT=43200
DT=1800
STARTTIME=0
ENDTIME=$BIGDT

# delete all old runs
rm -rf [0-9]*

# set up the initial conditions on the initial mesh
imesh=0
mesh=${meshes[$imesh]}

# create case
newCase=$imesh
mkdir $newCase
cp -r blank/* $newCase

# Set the start and end time
sed -i -e 's/STARTTIME/'$STARTTIME'/g' -e 's/ENDTIME/'$ENDTIME'/g' \
    -e 's/WRITEINT/'$BIGDT'/g' -e 's/DT/'$DT'/g' \
    $newCase/system/controlDict

# link to the initial mesh
ln -sf $baseMesh/$mesh/0/polyMesh $newCase/constant/polyMesh
ln -sf ../$baseMesh/$mesh/0/dualMesh/polyMesh $newCase/constant/dualMesh/polyMesh

# create mountain and initial conditions
setWilliMountain -case $newCase
setWilliMountain -case $newCase -region dualMesh
setWilliSteady -case $newCase

# initial run before mesh adaptivity
shallowWaterFoam -case $newCase |& tee $newCase/log

# calculate vorticity and divergence of the solution
vorticityDual -case $newCase

# map to new meshes and continue the simulation
imesh=1
while [ $imesh != ${#meshes[*]} ]; do
    case=$newCase
    newCase=$imesh
    mesh=${meshes[$imesh]}
    let STARTTIME=$STARTTIME+$BIGDT
    let ENDTIME=$ENDTIME+$BIGDT
    echo case $case newCase $newCase mesh $mesh imesh $imesh \
         from $STARTTIME to $ENDTIME
    
    mkdir $newCase
    cp -r blank/* $newCase
    
    # Set the start and end time
    sed -i -e 's/STARTTIME/'$STARTTIME'/g' -e 's/ENDTIME/'$ENDTIME'/g' \
        -e 's/WRITEINT/'$BIGDT'/g' -e 's/DT/'$DT'/g' \
        $newCase/system/controlDict

    # link to new mesh
    ln -sf $baseMesh/$mesh/0/polyMesh $newCase/constant/polyMesh
    ln -sf ../$baseMesh/$mesh/0/dualMesh/polyMesh $newCase/constant/dualMesh/polyMesh
    
    # create mountain and initial conditions
    setWilliMountain -case $newCase
    setWilliMountain -case $newCase -region dualMesh

    # map height, vorticity and divergence onto new mesh
    mapFields -consistent -case $newCase -fields '(h divhu)' -noLagrangian \
              -sourceTime latestTime $case
    mapFields -consistent -case $newCase -fields '(pv)' -noLagrangian \
        -sourceTime latestTime -sourceRegion dualMesh -targetRegion dualMesh \
        $case
     
    # invert vorticity and divergence on the new mesh to calculate velocity
    invertVorticityDivergence -case $newCase

    # debugging
    if [ $imesh == 1 ] || [ $imesh == 2 ]; then
        vorticityDual -case $newCase -time $STARTTIME
        plotPatchData -time $STARTTIME -case $newCase -region dualMesh vorticity
        gv $newCase/$STARTTIME/vorticity.eps &
        plotPatchData -time $STARTTIME -case $newCase -region dualMesh pv
        gv $newCase/$STARTTIME/pv.eps &
        plotPatchData -time $STARTTIME -case $newCase h
        gv $newCase/$STARTTIME/h.eps &
    fi
        
    # start the new simulation
    shallowWaterFoam -case $newCase |& tee $newCase/log
    
    # vorticity and divergence for the newCase
    vorticityDual -case $newCase

    let imesh=$imesh+1
done

### plots on new mesh
#time=$STARTTIME
#time=$ENDTIME
#plotPatchData -case $newCase -time $time h
#gv $newCase/$time/h.eps &
#plotPatchData -case $newCase -time $time divhu
#gv $newCase/$time/divhu.eps &
#plotPatchData -case $newCase -time $time -region dualMesh pv
#gv $newCase/$time/pv.eps &
#plotPatchData -case $newCase -time $time hU
#gv $newCase/$time/hU.eps &

