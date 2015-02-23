#!/bin/bash -ve

#gmtPlot runAll/plotAllDiags

if [ "$#" -eq 0 ]; then
    time=432000
else
    time=$1
fi

# h and U errors
#sumFields $time Udiff $time Uf 0 Uf -scale0 1 -scale1 -1
sumFields $time hDiff $time h 0 h  -scale0 1 -scale1 -1
plotPatchData -time $time hDiff15
gv $time/hDiff15.eps &

#plotPatchData hU -time $time
#ev $time/hU.eps

#sumFields -region dualMesh $time pvDiff $time pv 0 pv -scale0 1 -scale1 -1
#plotPatchData -time $time -region dualMesh pvDiff
#gv $time/pvDiff.eps &

#plotPatchData hU
#eps2gif hU.gif 0/hU.eps ????/hU.eps ?????/hU.eps ??????/hU.eps

