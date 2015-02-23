/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshWithDual, 0);
}

// * * * * * * * * * * * Constructor Helper funciton * * * * * * * * * * * * //
const Foam::polyPatch& Foam::fvMeshWithDual::findBottomPatch
(
    const word patchName
) const
{
    label patchID = boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        return boundaryMesh().last();
    }
    return boundaryMesh()[patchID];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshWithDual::fvMeshWithDual
(
    const IOobject& io,
    const bool readOverrides
)
:
    fvMesh(io),
    earthProperties_
    (
        IOobject
        (
            "earthProperties", time().constant(), *this, IOobject::READ_IF_PRESENT
        )
    ),
    isSpherical_
    (
        earthProperties_.found("earthRadius")
    ),
    patchName_(earthProperties_.lookupOrDefault<word>("surfacePatch", "originalPatch")),
    bottomPatch_(findBottomPatch(patchName_)),
    Omega_(earthProperties_.lookupOrDefault<dimensionedVector>("Omega", dimensionedVector("O", dimless, vector(0.,0.,0.)))),
    magg_(earthProperties_.lookupOrDefault<dimensionedScalar>("magg", dimensionedScalar("magg", dimless, scalar(0)))),
    OmegaHat_(unitVector(Omega_.value())),
    earthRadius_(earthProperties_.lookupOrDefault<dimensionedScalar>("earthRadius", dimensionedScalar("earthRadius", dimLength, scalar(1)))),
    nLevels_(nCells()/(max(bottomPatch_.size(),1))),
    faceToPatchEdgePtr_(NULL),
    faceToPatchFacePtr_(NULL),
    pointToPatchPointPtr_(NULL),
    cellToPatchFacePtr_(NULL),
    dualMesh_(NULL),
    dualFaceMapPtr_(NULL),
    signFaceMapPtr_(NULL),
    faceVolOwnPtr_(NULL),
    faceVolNeiPtr_(NULL),
    faceVolPtr_(NULL),
    depthfPtr_(NULL),
    fpeAreaFracPtr_(NULL),
    facePointAreaFracPtr_(NULL),
    facePointsPtr_(NULL),
    edgeAreaOwnPtr_(NULL),
    edgeAreaNeiPtr_(NULL),
    HdiagPtr_(NULL),
    intersectionsPtr_(NULL),
    cellEdgeVolsPtr_(NULL),
    cellEdgeCellsPtr_(NULL),
    rHatPtr_(NULL),
    lonHatPtr_(NULL),
    latHatPtr_(NULL),
    rHatfPtr_(NULL),
    lonHatfPtr_(NULL),
    latHatfPtr_(NULL),
    lonPtr_(NULL),
    latPtr_(NULL),
    heightPtr_(NULL),
    lonfPtr_(NULL),
    latfPtr_(NULL),
    heightfPtr_(NULL),
    idirPtr_(NULL),
    jdirPtr_(NULL),
    kdirPtr_(NULL),
    ddirPtr_(NULL)
{
    if (isSpherical_) Info << "Creating spherical mesh" << endl;
    else Info << "Creating Cartesian mesh with dual" << endl;
    newFaceCentresAndAreas(readOverrides);
    newCellCentresAndVols(readOverrides);
    polyMesh::boundary_.updateMesh();
    polyMesh::boundary_.calcGeometry();
}

Foam::fvMeshWithDual::fvMeshWithDual
(
    const IOobject& io,
    const fvMeshWithDual& dualMesh__,
    const bool readOverrides
)
:
    fvMesh(io),
    earthProperties_
    (
        IOobject
        (
            "earthProperties", time().constant(), *this, IOobject::READ_IF_PRESENT
        )
    ),
    isSpherical_
    (
        earthProperties_.found("earthRadius")
    ),
    patchName_(earthProperties_.lookupOrDefault<word>("surfacePatch", "originalPatch")),
    bottomPatch_(findBottomPatch(patchName_)),
    Omega_(earthProperties_.lookupOrDefault<dimensionedVector>("Omega", dimensionedVector("O", dimless, vector(0.,0.,0.)))),
    magg_(earthProperties_.lookupOrDefault<dimensionedScalar>("magg", dimensionedScalar("magg", dimless, scalar(0)))),
    OmegaHat_(unitVector(Omega_.value())),
    earthRadius_(earthProperties_.lookupOrDefault<dimensionedScalar>("earthRadius", dimensionedScalar("earthRadius", dimLength, scalar(1)))),
    nLevels_(nCells()/(max(bottomPatch_.size(),1))),
    faceToPatchEdgePtr_(NULL),
    faceToPatchFacePtr_(NULL),
    pointToPatchPointPtr_(NULL),
    cellToPatchFacePtr_(NULL),
    dualMesh_(&dualMesh__),
    dualFaceMapPtr_(NULL),
    signFaceMapPtr_(NULL),
    faceVolOwnPtr_(NULL),
    faceVolNeiPtr_(NULL),
    faceVolPtr_(NULL),
    depthfPtr_(NULL),
    fpeAreaFracPtr_(NULL),
    facePointAreaFracPtr_(NULL),
    facePointsPtr_(NULL),
    edgeAreaOwnPtr_(NULL),
    edgeAreaNeiPtr_(NULL),
    HdiagPtr_(NULL),
    intersectionsPtr_(NULL),
    cellEdgeVolsPtr_(NULL),
    cellEdgeCellsPtr_(NULL),
    rHatPtr_(NULL),
    lonHatPtr_(NULL),
    latHatPtr_(NULL),
    rHatfPtr_(NULL),
    lonHatfPtr_(NULL),
    latHatfPtr_(NULL),
    lonPtr_(NULL),
    latPtr_(NULL),
    heightPtr_(NULL),
    lonfPtr_(NULL),
    latfPtr_(NULL),
    heightfPtr_(NULL),
    idirPtr_(NULL),
    jdirPtr_(NULL),
    kdirPtr_(NULL),
    ddirPtr_(NULL)
{
    if (isSpherical_) Info << "Creating spherical mesh" << endl;
    else Info << "Creating Cartesian mesh with dual" << endl;
    newFaceCentresAndAreas(readOverrides);
    newCellCentresAndVols(readOverrides);
    polyMesh::boundary_.updateMesh();
    polyMesh::boundary_.calcGeometry();
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshWithDual::newFaceCentresAndAreas
(
    const bool readOverrides
)
{
    // check for faceCentres and faceAreas in mesh directory to overwrite
    IOobject newFaceCentresHeader
    (
        "faceCentres", time().findInstance(meshDir(), "points"),
        meshSubDir, *this, IOobject::MUST_READ, IOobject::NO_WRITE
    );
    IOobject newFaceAreasHeader
    (
        "faceAreas", time().findInstance(meshDir(), "points"),
        meshSubDir, *this, IOobject::MUST_READ, IOobject::NO_WRITE
    );

    // We can either read both, none or faceCentres only.
    
    // neither are present or not needed so calculate and write out if needed
    if
    (
        (!readOverrides
     || (!newFaceCentresHeader.headerOk() && !newFaceAreasHeader.headerOk()))
     && isSpherical()
    )
    {
        Pout<< "fvMeshWithDual::readNewFaceCentres() : "
            << "Calculating face centres and face areas" << endl;

        pointIOField fCtrs
        (
            IOobject
            (
                "faceCentres", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            //pointField(nFaces(), point::zero)
            faceCentres()
        );
        pointIOField fAreas
        (
            IOobject
            (
                "faceAreas", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            pointField(nFaces(), point::zero)
        );
        
        calcSphFaceCentresAndAreas(points(), fCtrs, fAreas);
        overrideFaceCentres(fCtrs);
        overrideFaceAreas(fAreas);
        if (readOverrides)
        {
            Info << "Writing face centres and areas" << endl;
            fCtrs.write();
            fAreas.write();
        }
    }
    // We can either read face areas but not face centres.
    else if
    (
        !newFaceCentresHeader.headerOk() && newFaceAreasHeader.headerOk()
    )
    {
        FatalErrorIn("fvMeshWithDual::readNewFaceCentresAndAreas")
            << " cannot read just \n"
            << newFaceAreasHeader.objectPath() << " without\n"
            << newFaceCentresHeader.objectPath()
            << exit(FatalError);
    }
    // we can read both
    else if (newFaceCentresHeader.headerOk() && newFaceAreasHeader.headerOk())
    {
        Pout << "Over-riding calculated face centres with\n"
             << newFaceCentresHeader.objectPath() << "\nand face areas with\n"
             << newFaceAreasHeader.objectPath() << endl;
        const pointIOField newFaceCentres(newFaceCentresHeader);
        const pointIOField newFaceAreas(newFaceAreasHeader);
        overrideFaceCentres(newFaceCentres);
        overrideFaceAreas(newFaceAreas);
    }
    // just face centres present
    else if (newFaceCentresHeader.headerOk())
    {
        Pout << "Over-riding calculated face centres with\n"
             << newFaceCentresHeader.objectPath() << endl;

        pointIOField newFaceCentres(newFaceCentresHeader);
        
        pointIOField fAreas
        (
            IOobject
            (
                "faceAreas", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            pointField(nFaces())
        );
        
        if (isSpherical())
        {
            calcSphFaceCentresAndAreas
            (
                points(), newFaceCentres, fAreas, false
            );
            overrideFaceAreas(fAreas);
            Info << "Writing face areas" << endl;
            fAreas.write();
            //newFaceCentres.write();
        }
        overrideFaceCentres(newFaceCentres);
    }
}


void Foam::fvMeshWithDual::newCellCentresAndVols
(
    const bool readOverrides
)
{
    // check for cellCentres and cellVolumes in mesh directory to overwrite
    IOobject newCellCentresHeader
    (
        "cellCentres", time().findInstance(meshDir(), "points"),
        meshSubDir, *this, IOobject::MUST_READ, IOobject::NO_WRITE
    );
    IOobject newCellVolumesHeader
    (
        "cellVolumes", time().findInstance(meshDir(), "points"),
        meshSubDir, *this, IOobject::MUST_READ, IOobject::NO_WRITE
    );
    
    // We can either read both, none or cellCentres only
    // neither are present or not needed so calculate and write out if needed
    if
    (
        (!readOverrides
     || (!newCellCentresHeader.headerOk()&&!newCellVolumesHeader.headerOk()))
     && isSpherical()
    )
    {
        Pout<< "fvMeshWithDual::newCellCentresAndVols() : "
            << "Calculating cell centres and volumes"<< endl;

        pointIOField cCtrs
        (
            IOobject
            (
                "cellCentres", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            pointField(nCells(), point::zero)
            //cellCentres()
        );
        scalarIOList vols
        (
            IOobject
            (
                "cellVolumes", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            scalarList(nCells(), scalar(0))
        );
        
        calcSphCellCentresAndVols(faceCentres(), faceAreas(),cCtrs,vols);
        overrideCellCentres(cCtrs);
        overrideCellVols(vols);
        if (readOverrides)
        {
            Info << "Writing cell centres and volumes" << endl;
            cCtrs.write();
            vols.write();
        }
    }
    // We cannot read cell volumes without cell centres
    else if (!newCellCentresHeader.headerOk()&&newCellVolumesHeader.headerOk())
    {
        FatalErrorIn("fvMeshWithDual::newCellCentresAndVols")
            << " cannot read just "
            << newCellVolumesHeader.objectPath() << " without\n"
            << newCellCentresHeader.objectPath()
            << exit(FatalError);
    }
    // reading both
    else if (newCellVolumesHeader.headerOk() && newCellCentresHeader.headerOk())
    {
        Pout << "Over-riding calculated cell centres with\n"
             << newCellCentresHeader.objectPath()
             << "\nand cell volumes with\n"
             << newCellVolumesHeader.objectPath() << endl;
        const pointIOField newCellCentres(newCellCentresHeader);
        const scalarIOList newCellVolumes(newCellVolumesHeader);
        overrideCellCentres(newCellCentres);
        overrideCellVols(newCellVolumes);
    }
    // just read cell centres
    else if (newCellCentresHeader.headerOk())
    {
        Pout << "Over-riding calculated cell centres with\n"
             << newCellCentresHeader.objectPath() << endl;
        pointIOField newCellCentres(newCellCentresHeader);
                
        scalarIOList vols
        (
            IOobject
            (
                "cellVolumes", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            scalarList(nCells(), scalar(0))
        );
        
        if (isSpherical())
        {
            calcSphCellCentresAndVols
            (
                faceCentres(), faceAreas(), newCellCentres, vols, false
            );
            overrideCellVols(vols);
            Info << "Writing cell volumes" << endl;
            vols.write();
            //newCellCentres.write();
        }

        overrideCellCentres(newCellCentres);
    }
}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshWithDual::setDual(const fvMeshWithDual& dm)
{
    dualMesh_ = &dm;
    
    // Check that the vertical face areas of the dual mesh to be consistent
    // with the primal
//    if (isSpherical())
//    {
//        pointIOField dualfAreas
//        (
//            IOobject
//            (
//                "faceAreas", dm.time().findInstance(dm.meshDir(), "points"),
//                dm.meshSubDir, dm
//            ),
//            dm.faceAreas()
//        );
//        
//        forAll(dualfAreas, faced)
//        {
//            label faceI = dm.dualFaceMap()[faced];
//            if (faceI != -1)
//            {
//                dualfAreas[faced] = unitVector(dualfAreas[faced])
//                                    *dm.depthf()[faced]/deltaCoeffs()[faceI];
//            }
//        }
//        
//        //dm.overrideFaceAreas(dualfAreas);
//        //dualfAreas.write();
//        
//        if (max(mag(dualfAreas - dm.faceAreas())) >SMALL*earthRadius().value())
//        {
//            FatalErrorIn("fvMeshWithDual::setDual")
//                << "primal and dual have incompatable face areas and "
//                << "deltaCoeffs by " << max(mag(dualfAreas - dm.faceAreas()))
//                << abort(FatalError);
//        }
//    }
}

// ************************************************************************* //
