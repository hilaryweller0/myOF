// read in the meshes
//Foam::fvMeshWithCentres mesh
//(
//    Foam::IOobject
//    (
//        Foam::fvMesh::defaultRegion,
//        runTime.timeName(),
//        runTime,
//        Foam::IOobject::MUST_READ
//    )
//);

Foam::fvMeshWithCentres dualMesh
(
    Foam::IOobject
    (
        "dualMesh",
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

// Identify the patches
//const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);
const polyPatch& dualPatch(dualMesh.boundaryMesh()["originalPatch"]);

// find map of faces from dual to primal mesh
IOList<label> primalDualFaceMap
(
    IOobject
    (
        "primalDualFaceMap",
        runTime.findInstance(dualMesh.meshDir(), "points"),
        dualMesh.meshDir(),
        runTime
    ),
    labelList(dualMesh.nInternalFaces(), -1)
);

for(label faceI = 0; faceI < dualMesh.nInternalFaces(); faceI++)
{
    const label i = dualMesh.faceOwner()[faceI];
    const label j = dualMesh.faceNeighbour()[faceI];

    const labelList& iEdges = ePatch.pointEdges()[i];
    const labelList& jEdges = ePatch.pointEdges()[j];
    
    // In the primal mesh, find the edge between points i and j
    for(label iei = 0; iei < iEdges.size() && primalDualFaceMap[faceI] == -1; iei++)
    {
        for(label jej = 0; jej < jEdges.size() &&primalDualFaceMap[faceI]==-1;jej++)
        {
            if (iEdges[iei] == jEdges[jej])
            {
                primalDualFaceMap[faceI] = iEdges[iei];
            }
        }
    }

    if(primalDualFaceMap[faceI]==-1 || primalDualFaceMap[faceI] >= dualMesh.nInternalFaces())
    {
        FatalErrorIn("polyDualPatch")
            << "Cannot find face in between primal mesh points "
            << i << " and " << j <<  " for dual face " << faceI << nl
            << "iEdges = " << iEdges << nl
            << "jEdges = " << jEdges << nl
            << exit(FatalError);
    }
}
primalDualFaceMap.write();

// find map of faces from primal to dual mesh
IOList<label> dualPrimalFaceMap
(
    IOobject
    (
        "dualPrimalFaceMap",
        runTime.findInstance(dualMesh.meshDir(), "points"),
        dualMesh.meshDir(),
        runTime
    ),
    labelList(dualMesh.nInternalFaces(), -1)
);

forAll(primalDualFaceMap, i)
{
    dualPrimalFaceMap[primalDualFaceMap[i]] = i;
}
dualPrimalFaceMap.write();

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

