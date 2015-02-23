/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    advectionFoam

Description
    Solves a transport equation for a passive scalar

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

volScalarField readT(Time& runTime, fvMesh& mesh)
{
    Info<< "Reading field T\n" << endl;

    return volScalarField
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}

volScalarField calculateDivUf(Time& runTime, fvMesh& mesh, surfaceScalarField& phi)
{
    Info<< "Calculating divergence of phi\n" << endl;

    return volScalarField
    (
        IOobject(
            "divUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi)
    );
}

surfaceScalarField readOrCalculatePhi(argList& args, Time& runTime, fvMesh& mesh)
{
    if (args.options().found("usePhi"))
    {
        Info<< "Reading field phi\n" << endl;

        return surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Reading field Uf\n" << endl;

        surfaceVectorField Uf
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Calculating face flux field phi\n" << endl;

        return surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Uf & mesh.Sf()
        );
    }
}

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("usePhi", "use 0/phi rather than calculating it from 0/Uf");
    Foam::argList::addBoolOption("leapfrog", "use leapfrog timestepping scheme rather than RK2");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField T = readT(runTime, mesh);
    surfaceScalarField phi = readOrCalculatePhi(args, runTime, mesh);
    volScalarField divUf = calculateDivUf(runTime, mesh, phi);
    divUf.write();

    Info<< "\nCalculating advection\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (args.options().found("leapfrog"))
        {
            T = T.oldTime().oldTime() - 2*runTime.deltaT()*fvc::div(phi,T);
        }
        else
        {
            for (int corr=0; corr < 3; corr++)
            {
                T = T.oldTime() - runTime.deltaT() * 0.5 *
                (
                    fvc::div(phi, T) + fvc::div(phi, T.oldTime())
                );
            }
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
