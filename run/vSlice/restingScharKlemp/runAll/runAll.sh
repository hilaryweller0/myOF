#!/bin/bash -e

# swap CPC for CFC
for dir in */save/*quadUpCPC; do
    newDir=`echo $dir | sed 's/CPC/CFC/g'`
    echo $dir $newDir
    mkdir $newDir
    cp -r $dir/0 $dir/constant $dir/system $newDir
    sed -i 's/CPC/CFC/g' $newDir/system/fvSchemes
done

for case in */save/*H_quadUpCFC; do
    echo $case
    #gedit -s $case/system/fvSchemes $case/system/fvSolution
    exnerFoamNonOrthogIMEXEX -case $case > $case/log 2>&1 &
done

for case in */save/*gradp_quadUpCFC; do
    echo $case
    #gedit -s $case/system/fvSchemes $case/system/fvSolution
    exnerFoamIMEXEX -case $case > $case/log 2>&1 &
done

