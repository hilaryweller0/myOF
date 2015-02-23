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
    oodles

Description
    Incompressible LES solver for flow in a channel.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/LESmodel/LESmodel.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "createAverages.H"
#   include "initContinuityErrs.H"
#   include "createGradP.H"

    const dimensionedScalar& dt = runTime.deltaT();

    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dt
    );

    volVectorField gradp = fvc::grad(p);
    volVectorField divB0 = sgsModel->divB(U) & U;
    volVectorField U0 = U;
    surfaceScalarField phi0 = phi;
    volVectorField divPhiU0 = fvc::div(phi, U, "div(phi,U)");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for(runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readMCSPPAControls.H"

#       include "CourantNo.H"

        sgsModel->correct();

        if (MCS)
        {
            // Multistep Advection using Adams-Bashforth

            dimensionedScalar dti = dt/nConvectionSteps;
            surfaceScalarField dphinp1 = phi - phi0;

            volVectorField divPhiU0 = 
                (scalar(nConvectionSteps - 1)/nConvectionSteps)
               *fvc::div(phi, U, "div(phi,U)")
              + (1.0/nConvectionSteps)
               *fvc::div(phi0, U0, "div(phi,U)");
            
            U0 = U;
            phi0 = phi;

            for (int i=0; i<nConvectionSteps; i++)
            {
                volVectorField divPhiU = fvc::div
                (
                    phi + (scalar(i)/nConvectionSteps)*dphinp1,
                    U,
                    "div(phi,U)"
                );

                U -= dti*(1.5*divPhiU - 0.5*divPhiU0);
                U.correctBoundaryConditions();
                divPhiU0 = divPhiU;
            }

            /*
            phi0 = phi;

            surfaceScalarField CoCoeff = 
                sqr(dti)*mesh.surfaceInterpolation::deltaCoeffs()/mesh.magSf();

            for (int i=0; i<nConvectionSteps; i++)
            {
                surfaceScalarField phii
                (
                    "phii",
                    phi + ((scalar(i) + 0.5)/nConvectionSteps)*dphinp1
                );

                surfaceVectorField Ucd = fvc::interpolate(U, "CD");
                surfaceVectorField Uud = fvc::interpolate(U, "UD");

                surfaceScalarField dtCof = CoCoeff*mag(phii);

                surfaceVectorField dtUif = (dti - dtCof)*Ucd + dtCof*Uud;

                U -= fvc::div(phii*dtUif);
                U.correctBoundaryConditions();
            }
            */
        }
        else
        {
            // Explicit Advection using Adams-Bashforth

            volVectorField divPhiU = fvc::div(phi, U);

            solve
            (
                fvm::ddt(U)
              + 1.5*divPhiU - 0.5*divPhiU0
            );

            divPhiU0 = divPhiU;
        }

        dimensionedScalar rdt = 1.0/dt;

        // Implicit stress using Crank-Nicolson
        {
            fvVectorMatrix divB = sgsModel->divB(U);

            solve
            (
                fvm::Sp(rdt, U) - rdt*U
              + 0.5*(divB + divB0)
              + gradp
              - flowDirection*gradP
            );

            divB0 = divB & U;
        }

        // Remove the pressure gradient contribution from U
        U += dt*gradp;

        // Flux trediction based on predicted velocity
        phi = (fvc::interpolate(U) & mesh.Sf())
            + fvc::ddtPhiCorr(rUA, U, phi);

        // Adjust flux to obet continuity
        adjustPhi(phi, U, p);

        // Pressure projection
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(p) == rdt*fvc::div(phi)
            );

            pEqn.setReference(pRefCell, pRefValue);

            if (nonOrth == nNonOrthCorr)
            {
                pEqn.solve(mesh.solver(p.name() + "Final"));
            }
            else
            {
                pEqn.solve(mesh.solver(p.name()));
            }

            if (nonOrth == nNonOrthCorr)
            {
                surfaceScalarField gradpf = pEqn.flux();
                phi -= dt*gradpf;

                //gradp = fvc::reconstruct(gradpf);
                gradp = (U - fvc::reconstruct(phi))/dt;

                U -= dt*gradp;
                U.correctBoundaryConditions();
            }
        }

#       include "continuityErrs.H"


        // Correct driving force for a constant mass flow rate

        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStar = 
            (flowDirection & U)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to 
        // adjust the average flow-rate to the correct value
        dimensionedScalar gragPplus = 
            (magUbar - magUbarStar)/rUA.weightedAverage(mesh.V());

        U += flowDirection*rUA*gragPplus;

        gradP += gragPplus;

        Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
            << "pressure gradient = " << gradP.value() << endl;

#       include "calculateAverages.H"

        runTime.write();

#       include "writeNaveragingSteps.H"

#       include "writeGradP.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
