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

\*---------------------------------------------------------------------------*/

#include "TRiSK3dreconstructVec.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK3d
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> reconstructVec(const surfaceScalarField& un)
{
    const fvMesh& mesh = un.mesh();

    tmp<surfaceVectorField> tV
    (
        new surfaceVectorField
        (
            IOobject(un.name()+"_3", un.instance(), un.mesh()),
            fvc::interpolate(fvc::reconstruct(un*mesh.magSf()))
        )
    );
    surfaceVectorField& V = tV();
    
    V += (un*mesh.magSf() - (V&mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

    return tV;
}


tmp<surfaceVectorField> reconstructVec(const tmp<surfaceScalarField>& tun)
{
    tmp<surfaceVectorField> tU
    (
        TRiSK3d::reconstructVec(tun())
    );
    tun.clear();
    return tU;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK3d

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
