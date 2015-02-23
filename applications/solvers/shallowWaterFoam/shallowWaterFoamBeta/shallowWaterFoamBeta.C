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
    shallowWaterFoamFluxBeta

Description
    Transient Solver for multi-layer SW equation in staggered flux form on the
    beta plane

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
//#include "OFstream.H"
#include "fvCFD.H"
#include "Hoperator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    Foam::Info
        << "Create mesh with dual for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        false
    );
    #include "readSolutionDict.H"
    #include "readEarthProperties.H"
    #define dt runTime.deltaT()
    #define Hdiag mesh.Hdiag()
    #include "createFields.H"
//    #include "initialInvariants.H"
//    #include "invariants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
#       include "CourantNo.H"

        for (int corr=0; corr<nCorr; corr++)
        {
            for(int icorr = 0; icorr < nNonOrthCorr; icorr++)
            {
                dvSdt = -mesh.magSf()*(mesh.ddir() & ((F ^ Uf)*hf));
                // Non-linear part of momentum change
                if (!linear)
                {
                    dvSdt -= mesh.magSf()*(mesh.ddir() & 
                    (
                       fvc::interpolate(fvc::div(phi, U))
                    ));
                }
            
                // Momentum change without d direction pressure gradient
                hvS = hvS.oldTime()
                    + dt*(1-offCentre)*dvSdt.oldTime() + dt*offCentre*dvSdt;
                phi = Hops::ddirToFlux(hvS)
                    - dt*offCentre*magg*hf*Hops::ddirToFluxOffDiag
                      (
                          mesh.magSf()*fvc::snGrad(h)
                      );

                fvScalarMatrix hEqn
                (
                    fvm::ddt(h)
                  + fvc::div((1-offCentre)*phi.oldTime() + offCentre*phi)
                  - fvm::laplacian(sqr(offCentre)*magg*hf*dt*Hdiag, h)
                );

                if (corr < nCorr-1) hEqn.solve();
                else hEqn.solve(mesh.solver(h.name() + "Final"));

                if (icorr == nNonOrthCorr-1)
                {
                    phi += 1./offCentre*hEqn.flux();
                }
            }

            hvS -= offCentre*dt*mesh.magSf()*hf*magg*fvc::snGrad(h);

            // update full velocity and non-linear terms
            if (!linear) hf = fvc::interpolate(h);
            
//            Uf = Hops::faceReconstruct(phi)/hf;
//            U = Hops::cellReconstruct(phi)/h;
            Uf = Hops::faceReconstruct(phi/hf);
            U = Hops::cellReconstruct(phi/hf);
            U.correctBoundaryConditions();
        }
        
        // update dvSdt for the next time-step
        dvSdt = -mesh.magSf()*
        (
            (mesh.ddir() & ((F ^ Uf)*hf))
          + magg*hf*fvc::snGrad(h)
        );
        if (!linear)
        {
            dvSdt -= mesh.magSf()*
            (
                (mesh.ddir() & (fvc::interpolate(fvc::div(phi, U))))
            );
        }

        runTime.write();
//        #include "invariants.H"

//        // Velocity based on flux
//        Uf.write();
//        Uf = fvc::interpolate(fvc::reconstruct(phi))/hf;
//        Uf.rename("Ufphi");
//        Uf.write();
//        FatalErrorIn("") << abort(FatalError);

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
