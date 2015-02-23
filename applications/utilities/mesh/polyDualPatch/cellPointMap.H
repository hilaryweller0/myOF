    // Map from primal cell centres to dual points
    IOList<label> cellToPointMap
    (
        IOobject
        (
            "cellToPointMap",
            runTime.findInstance(dualMesh.meshDir(), "points"),
            dualMesh.meshDir(),
            runTime
        ),
        labelList(mesh.nCells(),-1)
    );

    for(label ip = 0; ip < mesh.nCells(); ip++)
    {
        cellToPointMap[ip] = dualPatch.whichPoint(ip%dualPatch.nPoints());
    }
    cellToPointMap.write();

