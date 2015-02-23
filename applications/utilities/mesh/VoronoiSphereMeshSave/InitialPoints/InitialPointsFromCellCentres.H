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
    Foam::InitialPointsFromCellCentres

Description
    Read the initial points from an existing OpenFOAM case and provide
    interpolation onto new points

SourceFiles
    InitialPointsFromCellCentres.C

\*---------------------------------------------------------------------------*/

#ifndef InitialPointsFromCellCentres_H
#define InitialPointsFromCellCentres_H

#include "InitialPoints.H"
#include "meshWithDual.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class InitialPointsFromCellCentres Declaration
\*---------------------------------------------------------------------------*/

class InitialPointsFromCellCentres
:
    public InitialPoints
{
    // OpenFoam Case that the mesh is read from
    fileName case_;
    
    // Time database for the mesh
    Time runTime_;

    // OpenFoam mesh with cell and face centres overwritten
    fvMeshWithDual mesh_;
    
    // Flag to say if required resolution is globally uniform
    bool uniform_;
    
    // Store the uniform resolutution
    scalar uniformRes_;
    
    // required resolution at the cell centres (for non-uniform)
    volScalarField reqRes_;
    
public:

    //- Runtime type information
    TypeName("InitialPointsFromCellCentres");

    // Constructors

        //- Construct from components
        InitialPointsFromCellCentres
        (
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~InitialPointsFromCellCentres()
    {}


    // Member Functions

        // Access
        const pointField& points() const { return mesh_.cellCentres();}
        
        //- Helper function to interpolate a volScalarField onto an arbitrary point
        scalar interpolate(const volScalarField&, const point& p) const;

        //- Interpolate required resolution onto an arbitrary point
        scalar requiredResolution(const point& p) const
        {
            return uniform_ ? uniformRes_ : interpolate(reqRes_, p);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
