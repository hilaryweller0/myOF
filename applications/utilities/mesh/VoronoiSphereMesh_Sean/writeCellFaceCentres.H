// Read mesh extrusion properties from earthProperties
//IOdictionary earthProperties
//(
//    IOobject
//    (
//        "earthProperties",
//        runTime.constant(),
//        runTime,
//        IOobject::MUST_READ_IF_MODIFIED
//    )
//);

//extrudeModels::linearRadial radialExtrude(earthProperties);

// Writing the cell and face centres
const scalar r1 = radialExtrude.Rsurface();
const scalar r2 = radialExtrude.Router();
//const scalar Rmid = Foam::sqrt((sqr(r1) + r1*r2 + sqr(r2))/3.);
const scalar Rmid = 0.5*(r1 + r2);

//const polyPatch& dualPatch(dualMesh.boundaryMesh()["originalPatch"]);

// Declare the new cell centres
pointIOField cellCentres
(
    IOobject("cellCentres", runTime.timeName(), "polyMesh", mesh),
    cgalMesh.triPoints()
);

cellCentres *= Rmid/mag(cellCentres);
cellCentres.write();

// Declare the new dual cell centres
pointIOField dualCellCentres
(
    IOobject("cellCentres", runTime.timeName(), "polyMesh", dualMesh),
    pointField(dualMesh.nCells(), vector::zero)
);

forAll(dualCellCentres, celli)
{
    dualCellCentres[celli] = circumcentre
    (
        dualPatch.localPoints()[dualPatch[celli][0]],
        dualPatch.localPoints()[dualPatch[celli][1]],
        dualPatch.localPoints()[dualPatch[celli][2]]
    );
    dualCellCentres[celli] *= Rmid/mag(dualCellCentres[celli]);
}
dualCellCentres.write();

// find the cross-over points between Voronoi and Delaunay meshes
// and override the dual face centres
const edgeList& dualPatchEdges = dualPatch.edges();
pointIOField dualFaceCentres
(
    IOobject("faceCentres", runTime.timeName(), "polyMesh", dualMesh),
    pointField(dualPatchEdges.size() + 2*dualCellCentres.size())
);
forAll(dualPatchEdges, ie)
{
    // find the intersection between the points at either end of the
    // edge and the Delaunay edges
    dualFaceCentres[ie] = dualPatch.localPoints()[dualPatchEdges[ie][0]]
                        + dualPatch.localPoints()[dualPatchEdges[ie][1]];
    dualFaceCentres[ie] *= Rmid/mag(dualFaceCentres[ie]);
}

label m = dualPatchEdges.size();
forAll(dualCellCentres, i)
{
    dualFaceCentres[m + i] = r1/Rmid*dualCellCentres[i];
}
m = dualPatchEdges.size()+dualCellCentres.size();
forAll(dualCellCentres, i)
{
    dualFaceCentres[m + i] = r2/Rmid*dualCellCentres[i];
}

dualFaceCentres.write();

// Use the cross-over points for the primal face centres
pointIOField faceCentres
(
    IOobject("faceCentres", runTime.timeName(), "polyMesh", mesh),
    pointField(mesh.nFaces())
);
for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
{
    faceCentres[facei] = dualFaceCentres[dualPrimalFaceMap[facei]];
}

m = mesh.nInternalFaces();
forAll(cellCentres, i)
{
    faceCentres[m + i] = r1/Rmid*cellCentres[i];
}
m += cellCentres.size();
forAll(cellCentres, i)
{
    faceCentres[m + i] = r2/Rmid*cellCentres[i];
}

faceCentres.write();

Info << "Primal face 96 is at " << faceCentres[96] << nl
     << "dual face is at " << dualFaceCentres[dualPrimalFaceMap[96]] << endl;
Info << "Dual face 96 is at " << dualFaceCentres[96] << nl
     << "dual face is at " << faceCentres[primalDualFaceMap[96]] << endl;

