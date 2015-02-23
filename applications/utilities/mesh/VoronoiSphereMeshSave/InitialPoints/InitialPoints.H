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
    The Initial Voronoi generating points, from either a mesh or a list of points
    
Class
    Foam::InitialPoints

Description
    Abstract base class for points from a file or from a list of points.
    Also provides the required resolution and mesh density interpolated
    onto arbitrary points

SourceFiles
    initialPoints.C
    

\*---------------------------------------------------------------------------*/

#ifndef InitialPoints_H
#define InitialPoints_H

#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class InitialPoints Declaration
\*---------------------------------------------------------------------------*/

class InitialPoints
:
    public IOdictionary
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        InitialPoints(const InitialPoints&);

        //- Disallow default bitwise assignment
        void operator=(const InitialPoints&);


public:

    //- Runtime type information
    TypeName("InitialPoints");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            InitialPoints,
            dictionary,
            (
                const IOdictionary& dict
            ),
            (dict)
        );


    // Constructors
    
        //- Construct from dictionary
        InitialPoints(const IOdictionary& dict);
    
    // Selectors
    
        //- Return a reference to the selected InitialPoints
        static autoPtr<InitialPoints> New
        (
            const IOdictionary& dict
        );

    //- Destructor
    virtual ~InitialPoints(){};
    
    // Member Functions
    
        //- Return the list of points
        virtual const pointField& points() const = 0;

        //- Interpolate required resolution onto an arbitrary point
        virtual scalar requiredResolution(const point& p) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
