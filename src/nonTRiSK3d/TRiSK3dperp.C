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

#include "TRiSKperp.H"
#include "fvMesh.H"
#include "TRiSKData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceScalarField> perp(const surfaceScalarField& un)
{
    const fvMesh& mesh = un.mesh();

    const TRiSKData& triskData = TRiSKData::New(mesh);

    tmp<surfaceScalarField> treconField
    (
        new surfaceScalarField
        (
            IOobject(un.name()+"_perp", un.instance(), un.mesh()),
            un.mesh(),
            dimensionedScalar(un.name(), un.dimensions(), scalar(0))
        )
    );
    surfaceScalarField& v = treconField();
    
    forAll(v, faceI)
    {
        const labelList& stencil = triskData.stencil(faceI);
        const scalarList& weight = triskData.lbyd_weights(faceI);
        forAll(stencil, i)
        {
            v[faceI] -= un[stencil[i]]*weight[i];
        }
    }

    return treconField;
}


tmp<surfaceScalarField> perp(const tmp<surfaceScalarField>& tun)
{
    tmp<surfaceScalarField> tuperp
    (
        TRiSK::perp(tun())
    );
    tun.clear();
    return tuperp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
