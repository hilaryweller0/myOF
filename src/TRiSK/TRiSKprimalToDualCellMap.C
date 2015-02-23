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

#include "TRiSKprimalToDualCellMap.H"
#include "TRiSKconserveInterp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > primalToDualCellMap
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vf.mesh())));
    const fvMeshWithDual& dualMesh = mesh.dualMesh();

//    tmp<GeometricField<Type, fvPatchField, volMesh> > tvfd
//    (
//        new GeometricField<Type, fvPatchField, volMesh>
//        (
//            IOobject(vf.name(), vf.instance(), dualMesh),
//            faceToCellMap(mesh.dualMap(conserveInterp(vf)))
//        )
//    );
    
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvfd
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject(vf.name(), vf.instance(), dualMesh),
            dualMesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero),
            vf.boundaryField().types()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vfd = tvfd();
    
    forAll(vf, cellI)
    {
        const labelList& edgeCells = mesh.cellEdgeCells()[cellI];
        const scalarList& edgeVols = mesh.cellEdgeVols()[cellI];
        forAll(edgeCells, i)
        {
            label celld = edgeCells[i];
            vfd[celld] += edgeVols[i]*vf[cellI];
        }
    }
    vfd.internalField() /= dualMesh.V();
    
    return tvfd;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > primalToDualCellMap
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvfd
    (
        TRiSK::primalToDualCellMap(tvf())
    );
    tvf.clear();
    return tvfd;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
