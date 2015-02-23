/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Cell to face interpolation scheme. Included in fvMesh.

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "demandDrivenData.H"
#include "coupledFvPatch.H"
//Hilary Weller added for spherical geometry
#include "VectorSpaceFunctions.H"

void Foam::surfaceInterpolation::makeDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::makeDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    const fvMeshWithDual& meshwd 
         = *(dynamic_cast<const fvMeshWithDual*>(&(mesh_)));

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    deltaCoeffs_ = new surfaceScalarField
    (
        IOobject
        (
            "deltaCoeffs",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless/dimLength
    );
    surfaceScalarField& DeltaCoeffs = *deltaCoeffs_;

    if (meshwd.isSpherical())
    {
        DeltaCoeffs = 0.5*mesh_.magSf()/meshwd.faceVol();
    }
    else
    {
        // Set local references to mesh data
        const volVectorField& C = meshwd.C();
        const labelUList& owner = meshwd.owner();
        const labelUList& neighbour = meshwd.neighbour();

        forAll(owner, facei)
        {
            DeltaCoeffs[facei] = 1.0/mag(C[neighbour[facei]] - C[owner[facei]]);
        }

        forAll(DeltaCoeffs.boundaryField(), patchi)
        {
            DeltaCoeffs.boundaryField()[patchi] =
                1.0/mag(mesh_.boundary()[patchi].delta());
        }
    }
}


void Foam::surfaceInterpolation::makeNonOrthDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    nonOrthDeltaCoeffs_ = new surfaceScalarField
    (
        IOobject
        (
            "nonOrthDeltaCoeffs",
            mesh_.pointsInstance(),
            mesh_
        ),
        deltaCoeffs()
    );
}

// ************************************************************************* //
