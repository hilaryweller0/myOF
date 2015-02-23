

for dir in diamondCube/[1-9]*eq; do
    writeMeshObj -case $dir -patchFaces -time 0
    mv $dir/patchFaces_originalPatch_0.obj $dir/patch.obj
    cp $dir/patch.obj ${dir}_mesh.obj
    writeMeshObj -case $dir -patchFaces -time 0 -region dualMesh
    mv $dir/patchFaces_originalPatch_0.obj $dir/dualPatch.obj
    cp $dir/patch.obj ${dir}_dualMesh.obj
    rm $dir/patchFaces_otherSide_0.obj
done
