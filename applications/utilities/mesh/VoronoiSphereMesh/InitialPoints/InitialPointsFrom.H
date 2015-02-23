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
    Foam::InitialPointsFrom

Description
    Abstract base class for points from a file.
    Also provides the required resolution and mesh density interpolated
    onto arbitrary points

SourceFiles
    InitialPointsFrom.C
    

\*---------------------------------------------------------------------------*/

#ifndef InitialPointsFrom_H
#define InitialPointsFrom_H

#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "InitialPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class InitialPointsFrom Declaration
\*---------------------------------------------------------------------------*/

class InitialPointsFrom
:
    public IOdictionary,
    public InitialPoints
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        InitialPointsFrom(const InitialPointsFrom&);

        //- Disallow default bitwise assignment
        void operator=(const InitialPointsFrom&);


public:

    //- Runtime type information
    TypeName("InitialPointsFrom");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            InitialPointsFrom,
            dictionary,
            (
                const IOdictionary& dict
            ),
            (dict)
        );


    // Constructors
    
        //- Construct from dictionary
        InitialPointsFrom(const IOdictionary& dict);
    
    // Selectors
    
        //- Return a reference to the selected InitialPointsFrom
        static autoPtr<InitialPointsFrom> New
        (
            const IOdictionary& dict
        );

    //- Destructor
    virtual ~InitialPointsFrom(){};
    
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
