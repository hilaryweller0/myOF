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
    invertVorticityDivergence

Description
    Reads in height, pv and divhu on the primal grid and writes out Uf

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    word Htype
    (
        mesh.solutionDict().lookupOrDefault<word>("nonOrthogHtype", "diagonal")
    );

    const volScalarField f
    (
        IOobject("f", runTime.constant(), dualMesh),
        2*(dualMesh.Omega() & dualMesh.C())/mag(dualMesh.C())
    );
    
    Info << "mesh.Omega = " << mesh.Omega() << " dualMesh.Omega = " << dualMesh.Omega() << endl;
    
    const volScalarField fprimal
    (
        IOobject("f", runTime.constant(), mesh),
        2*(mesh.Omega() & mesh.C())/mag(mesh.C())
    );
    
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "Reading field pv on the dual mesh\n" << endl;
        volScalarField pv
        (
            IOobject("pv", runTime.timeName(), dualMesh, IOobject::MUST_READ),
            dualMesh
        );
        
        Info << "Reading field h\n" << endl;
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        
        Info << "Reading field divhu\n" << endl;
        volScalarField divhu
        (
            IOobject("divhu", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );

        Info << "Calculating pv, h and vorticity on the dual mesh\n" << endl;
        //volScalarField pv("pv",TRiSK::primalToDualCellMap(pvp));
        volScalarField hd("hd",TRiSK::primalToDualCellMap(h));
        volScalarField vorticity("vorticity", pv*hd - f);
        vorticity.write();
        

        Info << "Calculating the streamfunction and velocity potential\n" << endl;
        volScalarField streamFunc
        (
            IOobject("streamFunc", runTime.timeName(), dualMesh),
            dualMesh,
            dimensionedScalar("streamFunc", dimensionSet(0,2,-1,0,0), scalar(0))
        );
        volScalarField velPot
        (
            IOobject("velPot", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("velPot", dimensionSet(0,3,-1,0,0), scalar(0))
        );
        fvScalarMatrix streamFuncEqn(fvm::laplacian(streamFunc) == vorticity);
        streamFuncEqn.setReference(0, scalar(0));
        streamFuncEqn.solve();
        fvScalarMatrix velPotEqn(fvm::laplacian(velPot) == divhu);
        velPotEqn.setReference(0, scalar(0));
        velPotEqn.solve();

        Info << "Composing components of the velocity\n" << endl;
        surfaceScalarField hf = linearInterpolate(h);
        surfaceScalarField phi
             = dualMesh.dualFluxMap(fvc::snGrad(streamFunc)*dualMesh.magSf())*hf
             + fvc::snGrad(velPot)*mesh.magSf();

        surfaceVectorField Uf
        (
            "Uf",
            TRiSK::faceReconstruct(phi/(hf*mesh.magSf()), mesh.idir())
        );
        Uf.write();
    }
}


// ************************************************************************* //

