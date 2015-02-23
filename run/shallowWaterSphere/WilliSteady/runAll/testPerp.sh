#!/bin/bash -e

time=432000
for case in *C/*/save/* HRbucky/*/save/* ; do
    testPerp -case $case > $case/testPerpLog
    mv $case/0/Herror $case/$time
    #plotPatchData -case $case -time $time Herror
    #gv $case/$time/Herror.eps &
done

