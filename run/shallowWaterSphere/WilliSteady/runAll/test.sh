#!/bin/bash -e

testTRiSK

# The velocity on the primal mesh and mapped to the dual and back
# and other reconstructions
sumFields 0 Uerror 0 Urecon 0 Uf -scale0 1 -scale1 -1
plotPatchData Utest
gv 0/Utest.eps &


plotPatchData divu -region dualMesh
gv 0/divu.eps &
plotPatchData divuError -region dualMesh
gv 0/divuError.eps &
plotPatchData -region dualMesh curlu
gv 0/curlu.eps &
plotPatchData KE
gv 0/KE.eps &

