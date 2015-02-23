/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    exnerFoamIMEXEX

Description
    Transient Solver for buoyant, inviscid, incompressible, non-hydrostatic flow
    using a simultaneous solution of Exner, theta and phi
    Sub-stepping for advection (IMEXEX)

\*---------------------------------------------------------------------------*/

//#include "Hops.H"
#include "fvCFD.H"
#include "ExnerTheta.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
//    Hops H(mesh);
//    const surfaceScalarField gd("gd", g & H.delta());
    #define dt runTime.deltaT()
    #include "createFields.H"
    #include "initContinuityErrs.H"
    const dimensionedScalar initHeat = fvc::domainIntegrate(theta*rho);
    #include "initEnergy.H"
    #include "energy.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    const Switch SIgravityWaves(mesh.schemesDict().lookup("SIgravityWaves"));
    const int nAdvectionSubSteps = readLabel(itsDict.lookup("nAdvectionSubSteps"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "compressibleCourantNo.H"

        // update old time variables for Crank-Nicholson
        phi.oldTime() += (1-offCentre)*dt*dPhidt;

        // Advection sub stepping
        #include "advectionSubSteps.H"

        for (int ucorr=0; ucorr<nOuterCorr; ucorr++)
        {
            // update density according to the continuity equation
            solve
            (
                fvm::ddt(rho) + (1-offCentre)*divPhi + offCentre*fvc::div(phi)
            );
            // theta equation
            fvScalarMatrix thetaEqn
            (
                fvm::ddt(rho, theta)
              + (1-offCentre)*divPhitheta.oldTime()
              + offCentre*fvc::div(phi, theta)
            );
            thetaEqn.solve();

            // Exner and momentum equations
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "exnerEqn.H"
            }
        }
        
        // Final solutions for rho (and theta)
        solve(fvm::ddt(rho) + (1-offCentre)*divPhi + offCentre*fvc::div(phi));
        // theta equation
        fvScalarMatrix thetaEqn
        (
            fvm::ddt(rho, theta)
          + (1-offCentre)*divPhitheta.oldTime()
          + offCentre*fvc::div(phi, theta)
        );
        thetaEqn.solve();
        
        // Update rates of change for next time step
        thetaf = fvc::interpolate(theta);
        rhof = fvc::interpolate(rho);
        dPhidt += rhof*(gSf - mesh.magSf()*Cp*thetaf*fvc::snGrad(Exner)) - muSponge*phi;
        divPhi = fvc::div(phi);
        divPhitheta = fvc::div(phi, theta);
        
        #include "compressibleContinuityErrs.H"

        dimensionedScalar totalHeatDiff = fvc::domainIntegrate(theta*rho)
                                        - initHeat;
        Info << "Heat error = " << (totalHeatDiff/initHeat).value() << endl;
        #include "energy.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
