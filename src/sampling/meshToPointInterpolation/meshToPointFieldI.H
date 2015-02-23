/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class approxType>
template<class Type>
inline Field<Type> Foam::meshToPointField<approxType>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (stencils_.size() == 0)
    {
        FatalErrorIn("meshToPointField<approxType>::interpolate")
            << " points to interpolate onto not set"
            << exit(FatalError);
    }

    Field<Type> rVal(stencils_.size(), pTraits<Type>::zero);
    forAll(rVal, ip)
    {
        for(label is = 0; is < stencils_[ip].size(); is++)
        {
            rVal[ip] += weights_[ip][is]*vf[stencils_[ip][is]];
        }
    }
    
    return rVal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
