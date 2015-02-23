/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "TRiSK3dData.H"
#include "plane.H"
#include "sphericalGeometry.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::TRiSK3dData::patchName = "originalPatch";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::TRiSK3dData::calcAreasVols()
{
    const fvMesh& mesh = this->mesh();
    
    forAll(faceVolOwn_, faceI)
    {
        const label own = mesh.owner()[faceI];
        const label nei = mesh.neighbour()[faceI];
        faceVolOwn_[faceI] = mesh.Sf()[faceI]
                           & (mesh.Cf()[faceI] - mesh.C()[own]);
        faceVolNei_[faceI] = mesh.Sf()[faceI]
                           & (mesh.C()[nei] - mesh.Cf()[faceI]);
    }
    faceVolOwn_ *= 1./3.;
    faceVolNei_ *= 1./3.;
    faceVol_ = faceVolOwn_ + faceVolNei_;

    // Check that all volumes agree
    scalar faceVolSum = sum(faceVol_).value();
    scalar volSum     = sum(mesh.V()).value();
    if (mag(faceVolSum - volSum)/volSum > 1e-12)
    {
    Info //<< setprecision(12)
        << "Primal mesh point zero radius = " << mag(mesh.points()[0])
        <<"\n            last point radius = " << mag(mesh.points().last())
        <<"\n      cell centre zero radius = " << mag(mesh.C()[0])
        <<"\n     first face centre radius = " << mag(mesh.faceCentres()[0]) 
        <<"\n      last face centre radius = " << mag(mesh.faceCentres().last())
        << endl;

    WarningIn("TRiSK3dData::TRiSK3dData") //<< setprecision(12)
        << "faceVolSum = " << faceVolSum
        << " volSum = " << volSum
        << "\nnormalised diff (faceVolSum - volSum)/volSum = "
        << (faceVolSum - volSum)/volSum << endl;
//        << exit(FatalError);
    }
}

void Foam::TRiSK3dData::calcDeltaCoeffs()
{
    const fvMesh& mesh = this->mesh();

    // deltaCoeffs for the primal mesh
    surfaceScalarField deltaCoeffs
    (
        IOobject("deltaCoeffs", mesh.pointsInstance(), mesh.meshSubDir, mesh),
        mesh.magSf()/(2*faceVol_)
    );
    deltaCoeffs.write();
    
    Info << "Face  1496 deltaCoeffs = " << deltaCoeffs[1492] << " mesh.deltaCoeffs()[1492] = "
         << mesh.deltaCoeffs()[1492] << endl;
        
    // also need to do boundaries ...
}

void Foam::TRiSK3dData::checkDeltaCoeffs()
{
    const fvMesh& mesh = this->mesh();

    // Test deltaCoeffs
    scalarField deltaDiffNorm(mesh.nInternalFaces(), scalar(0));
    scalar maxDiff = 0;
    label maxAt = -1;
    forAll(deltaDiffNorm, faceI)
    {
        deltaDiffNorm[faceI] = 
        (mesh.deltaCoeffs()[faceI] - 0.5*mesh.magSf()[faceI]/faceVol_[faceI])
                                /mesh.deltaCoeffs()[faceI];
        if (mag(deltaDiffNorm[faceI]) > maxDiff)
        {
            maxDiff = mag(deltaDiffNorm[faceI]);
            maxAt = faceI;
        }
    }
    if (maxDiff > 2e-10)
    {
        FatalErrorIn("checkDeltaCoeffs")
            << "Maximum normalised difference between mesh_.deltaCoeffs() and mesh.magSf()/2V is " << maxDiff << " at face " << maxAt << " where mesh.deltaCoeffs() = "
            << mesh.deltaCoeffs()[maxAt] << " Sf = " << mesh.magSf()[maxAt]
            << " V = " << faceVol_[maxAt]
            << "\nS/2V = " << 0.5*mesh.magSf()[maxAt]/faceVol_[maxAt]
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TRiSK3dData, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TRiSK3dData::TRiSK3dData(const fvMesh& mesh)
:
    MeshObject<fvMesh, MoveableMeshObject, TRiSK3dData>(mesh),
    unitPatch_
    (
        mesh.boundaryMesh()[patchName].localFaces(),
        mesh.boundaryMesh()[patchName].localPoints()
        /mag(mesh.boundaryMesh()[patchName].localPoints())
    ),
    nLevels_(mesh.nCells()/unitPatch_.size()),
    rHat_(unitVector(mesh.C())),
    rHatf_(unitVector(mesh.Cf())),
    idir_(unitVector(mesh.Sf())),
    jdir_(unitVector(rHatf_ ^ idir_)),
    lonHat_(unitVector(vector(0,0,1) ^ rHat_)),
    latHat_(unitVector(rHat_ ^ lonHat_)),
    lonHatf_(unitVector(vector(0,0,1) ^ rHatf_)),
    latHatf_(unitVector(rHatf_ ^ lonHatf_)),
    faceVolOwn_
    (
        IOobject("faceVolOwn", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("V", dimVol, scalar(0))
    ),
    faceVolNei_
    (
        IOobject("faceVolNei", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("V", dimVol, scalar(0))
    ),
    faceVol_
    (
        IOobject("faceVol", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("V", dimVol, scalar(0))
    )
{
    if (debug)
    {
        Info<< "Constructing TRiSK3dData from fvMesh"
            << endl;
    }

    calcAreasVols();
    
//    // Check and/or calculate the delta coeffs
//    calcDeltaCoeffs();
//    checkDeltaCoeffs();
}


bool Foam::TRiSK3dData::movePoints()
{
    return false;
}

