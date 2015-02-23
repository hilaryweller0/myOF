    Info << "Creating extruded dual mesh" << endl;
    
    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    //const label nLayers = readLabel(earthProperties.lookup("nLayers"));
    
    extrudeModels::linearRadial radialExtrude(earthProperties);
    extrudedMesh dualMesh
    (
        IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime
        ),
        dualPatch,
        radialExtrude
    );
    
    // remove the empty zeroth boundary
    dualMesh.removeBoundary();
    DynamicList<polyPatch*> dualPatches(2);
    if (radialExtrude.nLayers() == 1)
    {
        dualPatches.append
        (
            new emptyPolyPatch
            (
                "originalPatch", dualPatch.size(),
                dualMesh.nInternalFaces(), 0, dualMesh.boundaryMesh(),
                "empty"
            )
        );
        dualPatches.append
        (
            new emptyPolyPatch
            (
                "otherSide", dualPatch.size(), 
                dualMesh.nInternalFaces() + dualPatch.size(),
                1,
                dualMesh.boundaryMesh(),
                "empty"
            )
        );
    }
    else
    {
        dualPatches.append
        (
            new polyPatch
            (
                "originalPatch", dualPatch.size(),
                dualMesh.nInternalFaces(), 0, dualMesh.boundaryMesh(),
                "patch"
            )
        );
        dualPatches.append
        (
            new polyPatch
            (
                "otherSide", dualPatch.size(), 
                dualMesh.nInternalFaces() + dualPatch.size(),
                1,
                dualMesh.boundaryMesh(),
                "patch"
            )
        );
    }

    dualMesh.addPatches(dualPatches);

