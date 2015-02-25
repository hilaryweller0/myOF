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
    departurePointStencil

Description
    Works out the stencil for the interpolation onto a departure point

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToPointField.H"
//#include "rbfFit.H"
#include "polyFit.H"
#include "centredCFCCellToCellStencilObject.H"
#include "centredCFCFCCellToCellStencilObject.H"
#include "centredCPCCellToCellStencilObject.H"

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
    #define dt runTime.deltaT()

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading field Uf\n" << endl;

    surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    surfaceScalarField psif
    (
        IOobject("psif", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("psif", dimless, scalar(0))
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

    volScalarField stencilWeights
    (
        IOobject("stencilWeights", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("stencilWeights", dimless, scalar(0))
    );

    // Departure points for each face
    surfaceVectorField departurePoints
    (
        IOobject("departurePoints", runTime.timeName(), mesh),
        mesh.Cf() - 2*Uf*dt
    );
    departurePoints.write();

    // Class for interpolating onto departure points
    const extendedCentredCellToCellStencil& meshStencil
        = centredCFCFCCellToCellStencilObject::New(mesh);
    meshToPointField<scalar, polyFit<oTWOPLUS>,  extendedCentredCellToCellStencil>
        meshToDep(mesh, departurePoints, meshStencil);
    
    // Loop through all cells and work out if each cell is a member of the 
    // stencil and what the weight is
    forAll(psi, cellI)
    {
        psi = 0;
        psi[cellI] = 1;
        psif.internalField() = meshToDep.interpolate(psi);
        if (psif[faceI] != 0)
        {
            stencil[cellI] = 1;
            stencilWeights[cellI] = psif[faceI];
        }
    }

    stencil.write();
    stencilWeights.write();
    psif = 0;
    psif[faceI] = 1;
    psif.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
