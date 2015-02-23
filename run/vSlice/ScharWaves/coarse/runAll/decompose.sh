#!/bin/bash -e

rm -rf processor*
decomposePar -constant
# correct processor boundaries of Uf
for file in processor*/0/Uf; do
    sed -i 's/(-10 0 0)/(10 0 0)/g' $file
done

# run in parallel and reconstruct
mpirun -np 4 exnerFoamH -parallel > log_p 2>&1
reconstructPar

