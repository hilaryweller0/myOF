# mesh generation

rm -rf [0-9]* constant/polyMesh
mkdir constant/polyMesh
cp constant/blockMeshDict constant/polyMesh
blockMesh

# pre processing
mkdir 0
cp init_0/T init_0/U 0
setFields

# run
rm -r 0.*
scalarTransportFoamExp

# post processing
for dir in 0.?; do
    rm -rf ${dir}0
    mv $dir ${dir}0
done
writeCellDataxyz T
paste 0/T.xyz 0.*/T.xyz | awk '{print $1, $4, $8, $16, $20, $24, $28, $32, $36, $40, $44, $48, $52}' > T.dat

#cp 0.10/T.xyz save/TlimitedCubic1.xyz
#cp 0.10/T.xyz save/TvanLeer.xyz
#cp 0.10/T.xyz save/Tcubic.xyz
#cp 0.10/T.xyz save/TcubicFit.xyz
#cp 0.10/T.xyz save/TcubicFitLimitedCubid.xyz

gmtPlot runAll/plotT.gmt
gmtinfo T.dat

#mv T.eps TvanLeer.eps
#mv T.eps TlimitedCubicVanLeer.eps
#mv T.eps TlimitedCubicUpwind.eps
#mv T.eps TcubicFit.eps
#mv T.eps Tupwind.eps
#mv T.eps TlimitedCubic0.eps
#mv T.eps Teps
#mv T.eps Tcubic.eps
#mv T.eps Tlinear.eps
#mv T.eps TlimitedLinear1.eps
#mv T.eps TcubicFitLimitedCubic.eps

#gmtPlot runAll/plotTcompare.gmt

