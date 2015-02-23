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

Class
    Foam::bellMountainEdge

Description
    bellMountainEdge class : defines a bell mountain shape in terms of a vector
    pointing to the top of the mountain and the radius

SourceFiles
    bellMountainEdge.C

\*---------------------------------------------------------------------------*/

#ifndef bellMountainEdge_H
#define bellMountainEdge_H

#include "curvedEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class bellMountainEdge Declaration
\*---------------------------------------------------------------------------*/

class bellMountainEdge
:
    public curvedEdge
{
    // Private data

        point p1_, p2_;
        vector H_;
        scalar a_;
        scalar A_;
        vector Hr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        bellMountainEdge(const bellMountainEdge&);

        //- Disallow default bitwise assignment
        void operator=(const bellMountainEdge&);

public:

    //- Runtime type information
    TypeName("bellMountain");


    // Constructors

        //- Construct from components
        bellMountainEdge
        (
            const pointField& points,
            const label start, const label end,
            const vector&, const scalar
        );

        //- Construct from Istream setting pointsList
        bellMountainEdge(const pointField& points, Istream&);


    // Destructor

        virtual ~bellMountainEdge(){}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
