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
    Reconstruct the perpendicular component of a vector field from the normal
    components at cell faces

SourceFiles
    TRiSKperp.C

\*---------------------------------------------------------------------------*/


#ifndef TRiSKperp_H
#define TRiSKperp_H

#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace TRiSK functions Declaration
\*---------------------------------------------------------------------------*/

namespace TRiSK
{
    tmp<surfaceScalarField> perp(const surfaceScalarField& uS);

    tmp<surfaceScalarField> perp(const tmp<surfaceScalarField>& tuS);

    tmp<surfaceVectorField> perpVec(const surfaceScalarField& uS);

    tmp<surfaceVectorField> perpVec(const tmp<surfaceScalarField>& tuS);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
