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

Namespace
    Foam

Description
    The Initial Voronoi generating points from a list of points
    
Class
    Foam::InitialPointsRaw

Description
    A pointField for initialising a VoronoiSphereMesh.
    Also provides the required resolution and mesh density interpolated
    onto arbitrary points

SourceFiles
    InitialPointsRaw.C
    

\*---------------------------------------------------------------------------*/

#ifndef InitialPointsRaw_H
#define InitialPointsRaw_H

#include "InitialPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class InitialPointsRaw Declaration
\*---------------------------------------------------------------------------*/

class InitialPointsRaw
:
    public InitialPoints
{
private:

    // Private data
    
        //- Hold reference to the pointField
        const pointField& points_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        InitialPointsRaw(const InitialPointsRaw&);

        //- Disallow default bitwise assignment
        void operator=(const InitialPointsRaw&);


public:

    //- Runtime type information
    TypeName("InitialPointsRaw");

    // Constructors
    
        //- Construct from dictionary
        InitialPointsRaw(const pointField& pts);
    
    //- Destructor
    virtual ~InitialPointsRaw(){};
    
    // Member Functions
    
        //- Return the list of points
        virtual const pointField& points() const { return points_;};

        //- Interpolate required resolution onto an arbitrary point
        virtual scalar requiredResolution(const point& p) const { return 1;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
