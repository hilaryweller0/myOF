/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    makeSK94perturb

Description
    Add a pertubation to theta following Skamarock and Klemp 1994

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    using namespace constant::mathematical;

    const volVectorField& C = mesh.C();

    Info << "\nReading initialProperties" << endl;

    IOdictionary initialProperties
    (
        IOobject
        (
            "initialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar thetaPrime(initialProperties.lookup("thetaPrime"));
    const dimensionedScalar bubbleCentre(initialProperties.lookup("bubbleCentre"));
    const dimensionedScalar halfWidth(initialProperties.lookup("bubbleHalfWidth"));
    const dimensionedScalar H(initialProperties.lookup("bubbleHeight"));

    Info<< "Reading theta_init\n" << endl;
    volScalarField theta_init
    (
        IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
        mesh
    );

    Info<< "Setting theta\n" << endl;
    volScalarField theta
    (
        IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
        theta_init
    );

    forAll(theta, celli)
    {
        theta[celli] += max
        (
            thetaPrime.value()*Foam::sin(pi*C[celli].z()/H.value())
            /(1 + sqr((C[celli].x() - bubbleCentre.value())
                /halfWidth.value())),
            scalar(0)
        );
    }

    Info<< "Writing theta\n" << endl;
    theta.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
