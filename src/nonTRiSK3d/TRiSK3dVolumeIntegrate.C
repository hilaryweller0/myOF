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

#include "TRiSKVolumeIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > volumeIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();
    const TRiSKData& triskData = TRiSKData::New(mesh);
    return triskData.faceVol()*sf.internalField();
}


template<class Type>
tmp<Field<Type> > volumeIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<Field<Type> > tf = TRiSK::volumeIntegrate(tsf());
    tsf.clear();
    return tf;
}


template<class Type>
dimensioned<Type>
domainIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + sf.name() + ')',
        dimVol*sf.dimensions(),
        gSum(TRiSK::volumeIntegrate(sf))
    );
}

template<class Type>
dimensioned<Type>
domainIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    dimensioned<Type> integral = domainIntegrate(tsf());
    tsf.clear();
    return integral;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
