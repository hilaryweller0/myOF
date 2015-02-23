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
    setScalarOverOrography

Description
    Set the velocity, U, and the initial scalar, T for scalar transport over
    orography

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

using namespace Foam::constant::mathematical;

scalar ScharCos(const scalar x, const scalar a) {
    return sqr(Foam::cos(M_PI*x/a));
}

scalar ScharCosSmooth(const scalar x, const scalar a, const scalar hm) {
    scalar h = 0;
    if (mag(x) < a)
    {
        h = hm*sqr(Foam::cos(0.5*M_PI*x/a));
    }
    return h;
}

Foam::scalar ScharExp(const scalar x, const scalar a, const scalar hm)
{
    return hm*Foam::exp(-sqr(x/a));
}

scalar height_schaerCos(scalar x, scalar a, scalar hm, scalar lambda) {
    return ScharCosSmooth(x, a, hm) * ScharCos(x, lambda);
}

scalar height_schaerExp(scalar x, scalar a, scalar hm, scalar lambda) {
    return ScharExp(x, a, hm) * ScharCos(x, lambda);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "x0", "int",
        "specify horizontal placement of tracer, overrides x0 specified in initialConditions dictionary"
     );
    Foam::argList::addOption
    (
        "tracerFieldFileName", "filename", 
        "specify the name of the tracer field file name (default 'T')"
    );
#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "Reading initial conditions" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "setScalarOverOrographyDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    IOdictionary velocityDict
    (
        IOobject
        (
            "velocityFieldDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Maximum wind speed
    const scalar u0(readScalar(velocityDict.lookup("maxVelocity")));
    // Initial maximum tracer value
    const scalar rho0(readScalar(initDict.lookup("rho0")));
    // Initial tracer position
    scalar x0(readScalar(initDict.lookup("x0")));
    const scalar z0(readScalar(initDict.lookup("z0")));
    // Half widths
    const scalar Ax(readScalar(initDict.lookup("Ax")));
    const scalar Az(readScalar(initDict.lookup("Az")));
    string tracerFieldFileName = "T";
    if (args.options().found("tracerFieldFileName"))
    {
        tracerFieldFileName = args.options()["tracerFieldFileName"];
    }

    if (args.options().found("x0"))
    {
        x0 = readScalar(IStringStream(args.options()["x0"])());
    }

    Info << "Creating initial tracer field " << tracerFieldFileName << endl;
    volScalarField T
    (
        IOobject(tracerFieldFileName, runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(tracerFieldFileName, dimless, scalar(0)),
        "zeroGradient"
    );
    
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        // Calculating T
        forAll(T, cellI)
        {
            const point& c = mesh.C()[cellI];
            
            // Centre of the tracer for this time step
            scalar x0t = x0 + u0*runTime.value();
        
            // Define r as used in the initial tracer field
            scalar r = Foam::sqrt(sqr((c.x()-x0t)/Ax)+sqr((c.z()-z0)/Az));
        
            if (r<=1)
            {
                T[cellI] = rho0*sqr(Foam::cos(M_PI*r/2));
            }
            else T[cellI] = 0;
        }
        T.correctBoundaryConditions();
        T.write();
    }
    
    Info<< "End\n" << endl;
    
    return(0);
}

// ************************************************************************* //
