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

#include "TRiSKconserveInterp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > conserveInterp
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vf.mesh())));
    
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject(vf.name(), vf.instance(), mesh),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
//            vf.boundaryField().types()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& fVo = mesh.faceVolOwn();
    const surfaceScalarField& fVn = mesh.faceVolNei();
    const surfaceScalarField& fV = mesh.faceVol();

    forAll(sf, faceI)
    {
        sf[faceI] = 
        (
            fVo[faceI]*vf[owner[faceI]]
          + fVn[faceI]*vf[neighbour[faceI]]
        )/fV[faceI];
    }
    
    forAll(mesh.boundary(), patchI)
    {
        sf.boundaryField()[patchI] = vf.boundaryField()[patchI];
    }

    return tsf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > conserveInterp
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        TRiSK::conserveInterp(tvf())
    );
    tvf.clear();
    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
