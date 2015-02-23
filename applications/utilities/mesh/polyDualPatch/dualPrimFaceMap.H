// find map of faces from dual to primal mesh
IOList<label> dualPrimalFaceMap
(
    IOobject
    (
        "dualFaceMap",
        runTime.findInstance(dualMesh.meshDir(), "points"),
        dualMesh.meshDir(),
        runTime
    ),
    labelList(dualMesh.nFaces(), -1)
);
const label nEdges = dualPatch.nEdges();
const label nCols = dualPatch.size();

// Create the map of dual to primal edges
const labelList edgeMap = findEdgeMap(dualPatch, ePatch);

// map edges of patch onto vertical faces
for(label ie = 0; ie < nEdges; ie++)
{
    for(label k = 0; k < radialExtrude.nLayers(); k++)
    {
        dualPrimalFaceMap[ie + k*(nEdges + nCols)]
             = edgeMap[ie] + k*(nEdges + ePatch.size());
    }
}

dualPrimalFaceMap.write();

// find map of faces from primal to dual mesh
IOList<label> primalDualFaceMap
(
    IOobject
    (
        "dualFaceMap",
        runTime.findInstance(mesh.meshDir(), "points"),
        mesh.meshDir(),
        runTime
    ),
    labelList(mesh.nFaces(), -1)
);

forAll(dualPrimalFaceMap, i)
{
    if (dualPrimalFaceMap[i] != -1)
    {
        primalDualFaceMap[dualPrimalFaceMap[i]] = i;
    }
}
primalDualFaceMap.write();

