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

#include "TRiSK3dke.H"
#include "TRiSK3dData.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK3d
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> ke(const surfaceScalarField& un)
{
    const fvMesh& mesh = un.mesh();

    tmp<volScalarField> tke
    (
        new volScalarField
        (
            IOobject("ke("+un.name()+")", un.instance(), mesh),
            mesh,
            dimensionedScalar("ke", sqr(un.dimensions()), scalar(0)),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& ke = tke();
    
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    forAll(own, facei)
    {
        scalar unSqr = sqr(un[facei]);
        const vector& Sf = mesh.Sf()[facei];
        const vector& Cf = mesh.Cf()[facei];
        ke[own[facei]] += unSqr*(Sf & (Cf - mesh.C()[own[facei]]));
        ke[nei[facei]] += unSqr*(Sf & (mesh.C()[nei[facei]] - Cf));
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<scalar>& pun = un.boundaryField()[patchi];
        const fvsPatchField<vector>& pSf = mesh.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& pCf = mesh.Cf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            const vector& Sf = pSf[facei];
            const vector& Cf = pCf[facei];
            ke[pFaceCells[facei]] += sqr(pun[facei])
                                 *(Sf & (Cf - mesh.C()[pFaceCells[facei]]));
        }
    }

    ke.internalField() /= 3*mesh.V();
    ke.correctBoundaryConditions();

    return tke;
}


tmp<volScalarField> ke(const tmp<surfaceScalarField>& tun)
{
    tmp<volScalarField> tke
    (
        TRiSK3d::ke(tun())
    );
    tun.clear();
    return tke;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK3d

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
