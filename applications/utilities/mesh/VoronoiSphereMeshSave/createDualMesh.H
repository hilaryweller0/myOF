// Read in the primal mesh and create the dual

Foam::fvMeshWithCentres mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

// Identify the patch to dualise and then extrude
const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);

// Create the dual points and faces
pointField dualPoints(ePatch.faceCentres());
faceList   dualFaces(ePatch.nPoints());
const labelListList& pointFaces = ePatch.pointFaces();
    
forAll(dualFaces, ip)
{
    const labelList& f = pointFaces[ip];
    dualFaces[ip].setSize(f.size());
    forAll(f, i)
    {
        dualFaces[ip][i] = f[i];
    }
    
    // Change order if necessary
    labelList& df = dualFaces[ip];
    for(bool noSwaps = false; !noSwaps;)
    {
        const point& a = dualPoints[df[0]];
        noSwaps = true;
        for(label i = 1; i < df.size()-1; i++)
        {
            const point& b = dualPoints[df[i]];
            const point& c = dualPoints[df[i+1]];            
            // Flip direction if necessary
            if ((((b-a)^(c-b)) & (a+b+c)) < 0) 
            {
                noSwaps = false;
                label dfi = df[i];
                df[i] = df[i+1];
                df[i+1] = dfi;
            }
        }
    }
}

PrimitivePatch<face, List, pointField> dualPPatch(dualFaces, dualPoints);

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
extrudedMesh dualEMesh
(
    IOobject
    (
        "dualMesh",
        runTime.timeName(),
        runTime
    ),
    dualPPatch,
    radialExtrude
);

// Give the dual mesh two empty boundaries
dualEMesh.removeBoundary();
DynamicList<polyPatch*> newPatches(2);
newPatches.append
(
//    new polyPatch
    new emptyPolyPatch
    (
        "originalPatch", dualPPatch.size(),
        dualEMesh.nInternalFaces(), 0, dualEMesh.boundaryMesh()
    )
);
newPatches.append
(
//    new polyPatch
    new emptyPolyPatch
    (
        "otherSide", dualPPatch.size(), 
        dualEMesh.nInternalFaces() + dualPPatch.size(), 1, dualEMesh.boundaryMesh()
    )
);
dualEMesh.addPatches(newPatches);

dualEMesh.write();

