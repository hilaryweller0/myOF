/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

InNamespace
    Foam::TRiSK

Description
    Calculate the curl on the dual mesh from the normal face velocities

SourceFiles
    TRiSKcurl.C

\*---------------------------------------------------------------------------*/


#ifndef TRiSKcurl_H
#define TRiSKcurl_H

#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace TRiSK functions Declaration
\*---------------------------------------------------------------------------*/

namespace TRiSK
{
    tmp<volScalarField> curl(const surfaceScalarField& vS);

    tmp<volScalarField> curl(const tmp<surfaceScalarField>& vS);

    tmp<surfaceVectorField> curl3d(const surfaceScalarField& uS);

    tmp<surfaceVectorField> curl3d(const tmp<surfaceScalarField>& uS);

    tmp<volScalarField> primalCurl(const surfaceScalarField& phi);

    tmp<volScalarField> primalCurl(const tmp<surfaceScalarField>& phi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
