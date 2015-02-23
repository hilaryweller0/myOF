/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    bellMountainEdge class :  defines a bell mountain shape in terms of a vector
    pointing to the top of the mountain and the radius

\*---------------------------------------------------------------------------*/

#include "bellMountainEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellMountainEdge, 0);
    addToRunTimeSelectionTable(curvedEdge, bellMountainEdge, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::bellMountainEdge::bellMountainEdge
(
    const pointField& points,
    const label start,
    const label end,
    const vector& H,
    const scalar a
)
:
    curvedEdge(points, start, end),
    p1_(points_[start_]),
    p2_(points_[end_]),
    H_(H),
    a_(a),
    A_(mag(p1_ - p2_)/a_),
    Hr_(p1_ - H_/(1 + sqr(0.5*A_)))
{}


// Construct from Istream
Foam::bellMountainEdge::bellMountainEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    p1_(points_[start_]),
    p2_(points_[end_]),
    H_(is),
    a_(readScalar(is)),
    A_(mag(p1_ - p2_)/a_),
    Hr_(p1_ - H_/(1 + sqr(0.5*A_)))
    //Hr_(p1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::bellMountainEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("bellMountainEdge::position(const scalar lambda) const")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    if (lambda < SMALL)
    {
        return p1_;
    }
    else if (lambda > 1-SMALL)
    {
        return p2_;
    }
    else
    {
        return Hr_ + lambda*(p2_ - p1_) + H_/(1 + sqr(A_*(lambda - 0.5)));
    }
}


Foam::scalar Foam::bellMountainEdge::length() const
{
    return mag(p1_ - p2_);
}


// ************************************************************************* //
