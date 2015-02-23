/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    calcPV

Description
    Calculates potentiail vorticity on the dual mesh from the vorticity
    on the dual and the height on the primal mesh

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    TRiSKData trisk(mesh);

    const volScalarField f
    (
        IOobject("f", runTime.constant(), dualMesh),
        2*(mesh.Omega() & dualMesh.C())/mag(dualMesh.C())
    );
    const volScalarField fprimal
    (
        IOobject("f", runTime.constant(), mesh),
        2*(mesh.Omega() & mesh.C())/mag(mesh.C())
    );
    
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "Reading field vorticity\n" << endl;
        volScalarField vorticity
        (
            IOobject("vorticity", runTime.timeName(), dualMesh, IOobject::MUST_READ),
            dualMesh
        );
        
        Info << "Reading field h\n" << endl;
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        
        volScalarField pv
        (
            "pv",
            (vorticity+f)/TRiSK::primalToDualCellMap(h)
        );
        pv.write();
    }
}


// ************************************************************************* //

