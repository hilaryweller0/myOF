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

#include "TRiSK3dinterpolateV.H"
#include "TRiSK3dData.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK3d
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> interpolateV(const volVectorField& U)
{
    const fvMesh& mesh = U.mesh();

    // References to the TRiSK3d data
    const TRiSK3dData& triskData = TRiSK3dData::New(mesh);

    volScalarField w(U & triskData.rHat());
//    volScalarField u(U & triskData.lonHat());
//    volScalarField v(U & triskData.latHat());
    volVectorField Uh(U - w*triskData.rHat());
    
    tmp<surfaceVectorField> tinterpolateV
    (
        new surfaceVectorField
        (
            IOobject("interpolateV("+U.name()+")", U.instance(), mesh),
            fvc::interpolate(Uh)
//            triskData.lonHatf()*fvc::interpolate(u)
//          + triskData.latHatf()*fvc::interpolate(v)
//          + triskData.rHatf()*fvc::interpolate(w)
        )
    );
    surfaceVectorField& Uf = tinterpolateV();
    
    Uf -= triskData.rHatf()*(Uf & triskData.rHatf());
    Uf += triskData.rHatf()*fvc::interpolate(w);
    
    return tinterpolateV;
}


tmp<surfaceVectorField> interpolateV(const tmp<volVectorField>& tU)
{
    tmp<surfaceVectorField> tinterpolateV
    (
        TRiSK3d::interpolateV(tU())
    );
    tU.clear();
    return tinterpolateV;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK3d

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
