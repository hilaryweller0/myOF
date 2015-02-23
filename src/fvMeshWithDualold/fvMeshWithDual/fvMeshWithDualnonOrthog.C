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
#include "fvc.H"
#include "midPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMeshWithDual::makeddir() const
{
    if (ddirPtr_) FatalErrorIn("fvMeshWithDual::makeddir")
                    << "ddir already exists" << abort(FatalError);
    ddirPtr_ = new surfaceVectorField("ddir", unitVector(delta()));
}


void fvMeshWithDual::makeHdiag() const
{
    if (HdiagPtr_) FatalErrorIn("fvMeshWithDual::makeHdiag")
                    << "Hdiag already exists" << abort(FatalError);

    HdiagPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "Hdiag", time().findInstance(meshDir(), "points"),
            meshSubDir, *this
        ),
        //dualMesh().signedDualMap(dualMesh().jdir()) & idir()
        ddir() & idir()
    );
    Hdiag().write();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
