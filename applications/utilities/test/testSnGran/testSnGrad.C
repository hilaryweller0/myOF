/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

Application
    testSnGrad

Description
    Calculates the snGrad of a given field

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    Foam::argList::validArgs.append("field");

    #include "setRootCase.H"
    #include "createTime.H"
    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createMesh.H"
    #include "orthogonalBoundaries.H"
    const word fieldName(args.args()[1]);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField vf
    (
        IOobject(fieldName, runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    surfaceScalarField snGradvf(word("snGrad"+fieldName), fvc::snGrad(vf));

    snGradvf.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
