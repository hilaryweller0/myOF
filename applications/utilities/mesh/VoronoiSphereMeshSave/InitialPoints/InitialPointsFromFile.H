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
    Foam::InitialPointsFromFile

Description
    Read the initial points from a file and provide constant interpolation
    onto points

SourceFiles
    initialPointsFromFile.C

\*---------------------------------------------------------------------------*/

#ifndef InitialPointsFromFile_H
#define InitialPointsFromFile_H

#include "InitialPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class InitialPointsFromFile Declaration
\*---------------------------------------------------------------------------*/

class InitialPointsFromFile
:
    public InitialPoints
{
    // points read in from file
    pointField points_;
    
    // globally uniform required resolution, and mesh density
    scalar reqRes_;
    
public:

    //- Runtime type information
    TypeName("InitialPointsFromFile");

    // Constructors

        //- Construct from components
        InitialPointsFromFile
        (
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~InitialPointsFromFile()
    {}


    // Member Functions

        // Access
        const pointField& points() const { return points_;}

        //- Required resolution at an arbitrary point
        scalar requiredResolution(const point& p) const { return reqRes_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
