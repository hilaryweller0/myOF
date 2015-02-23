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

#include "TRiSKFaceToCellMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > faceToCellMap
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(sf.mesh())));
    
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject(sf.name(), sf.instance(), mesh),
            mesh,
            dimensioned<Type>(sf.name(), sf.dimensions(), pTraits<Type>::zero),
            sf.boundaryField().types()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& fVo = mesh.faceVolOwn();
    const surfaceScalarField& fVn = mesh.faceVolNei();

    forAll(sf, faceI)
    {
        vf[owner[faceI]] += sf[faceI]*fVo[faceI];
        vf[neighbour[faceI]] += sf[faceI]*fVn[faceI];
    }
    
//    forAll(mesh.boundary(), patchI)
//    {
//        const labelUList& pFaceCells =
//            mesh.boundary()[patchI].faceCells();

//        Field<Type>& pvf = vf.boundaryField()[patchI];
//        const Field<Type>& psf = sf.boundaryField()[patchI];
//        const scalarField& pfVo = fVo.boundaryField()[patchI];

//        forAll(mesh.boundary()[patchI], facei)
//        {
//            pvf[pFaceCells[facei]] += psf[facei]*pfVo[facei];
//        }
//        
//        vf.boundaryField()[patchI] = sf.boundaryField()[patchI];
//    }

    vf.internalField() /= mesh.V();
    
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > faceToCellMap
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        TRiSK::faceToCellMap(tsf())
    );
    tsf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > faceToCellFluxMap
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(sf.mesh())));
    
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject(sf.name(), sf.instance(), mesh),
            mesh,
            dimensioned<Type>(sf.name(), sf.dimensions(), pTraits<Type>::zero),
            sf.boundaryField().types()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& fVo = mesh.faceVolOwn();
    const surfaceScalarField& fVn = mesh.faceVolNei();

    forAll(sf, faceI)
    {
        vf[owner[faceI]] += sf[faceI]*fVo[faceI];
        vf[neighbour[faceI]] -= sf[faceI]*fVn[faceI];
    }
    vf.internalField() /= mesh.V();
    
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > faceToCellFluxMap
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        TRiSK::faceToCellFluxMap(tsf())
    );
    tsf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
