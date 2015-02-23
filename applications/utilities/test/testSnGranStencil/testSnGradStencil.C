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
    getStencil

Description
    Works out the stencil for the given snGrad scheme for face faceI

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("faceI");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    const label faceI = readLabel(IStringStream(args.args()[1])());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    surfaceScalarField psif
    (
        IOobject("psif", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("psif", dimensionSet(0,-1,0,0,0), scalar(0))
    );

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("psi", dimless, scalar(0))
    );

    volScalarField stencil
    (
        IOobject("stencil", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("stencil", dimless, scalar(0))
    );

    // Loop through all cells and work out if each cell is a member of the 
    // stencil
    forAll(psi, cellI)
    {
        psi = 0;
        psi[cellI] = 1;
        psif = fvc::snGrad(psi);
        if (psif[faceI] != 0)
        {
            stencil[cellI] = 1;
        }
    }

    stencil.write();
    psif *= 0;
    psif[faceI] = 1;
    psif.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
