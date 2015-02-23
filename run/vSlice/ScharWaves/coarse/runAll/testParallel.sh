#!/bin/bash -e

rm -r [1-9]* processor*/[1-9]*
#mpirun -np 2 testParallelReconstructUf -parallel
mpirun -np 4 exnerFoamH -parallel
reconstructPar
mv 40 41
#testParallelReconstructUf
exnerFoamH

rm -rf 41/uniform
fields=`ls 41`
for field in $fields ; do
    sumFields 40 ${field}Diff 40 $field 41 $field -scale0 1 -scale1 -1
    read -p "Press [Enter] to continue..."
done

Writing 40/"gradPcoeffDiff"
Max = max(gradPcoeffDiff) [1 0 -1 0 0 0 0] 1.12057e-05 min = min(gradPcoeffDiff) [1 0 -1 0 0 0 0] -1.0252e-05
pyf         
Writing 40/"UDiff"
Max = max(UDiff) [1 0 -1 0 0 0 0] 5.21541e-07 min = min(UDiff) [1 0 -1 0 0 0 0] -5.96046e-07

plotPatchData -time 40 rhoDiff
gv 40/rhoDiff.eps &
plotPatchData -time 40 PsiDiff
gv 40/PsiDiff.eps &

plotPatchData -time 40 U
ev 40/U.eps
plotPatchData -time 40 thetaDiff
gv 40/thetaDiff.eps &
plotPatchData -time 40 ExnerDiff
gv 40/ExnerDiff.eps &

