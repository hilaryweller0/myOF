/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "extendedCentredCellToFaceStencil.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class TransformOp>
void Foam::extendedCentredCellToFaceStencil::collectData
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type> >& stencilFld,
    const TransformOp& top
) const
{
    // 1. Construct flat field for internal and patch data
    List<Type> flatFld(map().constructSize(), pTraits<Type>::zero);

    // Insert my internal values
    forAll(fld, cellI)
    {
        flatFld[cellI] = fld[cellI];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];

        label nCompact = pfld.patch().start();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    map().distribute(fld.mesh().globalData().globalTransforms(), flatFld, top);


    // 2. Pull to stencil
    stencilFld.setSize(elements_.size());

    forAll(elements_, elemI)
    {
        const labelList& cCells = elements_[elemI];
        const labelList& trafoCCells = transformedElements_[elemI];

        stencilFld[elemI].setSize(cCells.size()+trafoCCells.size());

        label nCompact = 0;
        forAll(cCells, i)
        {
            stencilFld[elemI][nCompact++] = flatFld[cCells[i]];
        }
        forAll(trafoCCells, i)
        {
            stencilFld[elemI][nCompact++] = flatFld[trafoCCells[i]];
        }
    }
}


template<class Type>
void Foam::extendedCentredCellToFaceStencil::collectData
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type> >& stencilFld
) const
{
    collectData(fld, stencilFld, mapDistribute::transform());
}


template<class Type>
void Foam::extendedCentredCellToFaceStencil::collectPositions
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type> >& stencilFld
) const
{
    collectData(fld, stencilFld, mapDistribute::transformPosition());
}


// ************************************************************************* //
