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

\*---------------------------------------------------------------------------*/

#include "InitialPointsFromPointField.H"
#include "addToRunTimeSelectionTable.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(InitialPointsFromPointField, 0);
addToRunTimeSelectionTable(InitialPoints, InitialPointsFromPointField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

InitialPointsFromPointField::InitialPointsFromPointField
(
    const IOdictionary& dict
)
:
    InitialPoints(dict),
    points_(0),
    reqRes_(dict.lookupOrDefault<scalar>("uniformResolution", scalar(1)))
{
    fileName ptFile(dict.lookup("initialCellCentresFile"));
    IOobject fileIO
    (
        ptFile, dict.db(), IOobject::MUST_READ
    );
    pointIOField newPoints(fileIO);
    points_ = pointField(Xfer<List<point> >(newPoints));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
