    // remove the empty zeroth boundary
    newMesh.removeBoundary();
    DynamicList<polyPatch*> newPatches(2);
    if (radialExtrude.nLayers() == 1)
    {
        newPatches.append
        (
            new emptyPolyPatch
            (
                "originalPatch", newPatch.size(),
                newMesh.nInternalFaces(), 0, newMesh.boundaryMesh(),
                "empty"
            )
        );
        newPatches.append
        (
            new emptyPolyPatch
            (
                "otherSide", newPatch.size(), 
                newMesh.nInternalFaces() + newPatch.size(),
                1,
                newMesh.boundaryMesh(),
                "empty"
            )
        );
    }
    else
    {
        newPatches.append
        (
            new polyPatch
            (
                "originalPatch", newPatch.size(),
                newMesh.nInternalFaces(), 0, newMesh.boundaryMesh(),
                "patch"
            )
        );
        newPatches.append
        (
            new polyPatch
            (
                "otherSide", newPatch.size(), 
                newMesh.nInternalFaces() + newPatch.size(),
                1,
                newMesh.boundaryMesh(),
                "patch"
            )
        );
    }

    newMesh.addPatches(newPatches);

