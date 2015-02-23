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
    icoSphereFoam

Description
    Transient Solver for incompressible Euler equations in 3d on the sphere
    with (non-)orthogonal TRiSK

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK.H"
//#include "OFstream.H"
#include "fvCFD.H"
#include "Hoperator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    TRiSKData trisk(mesh);
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

        // Momentum predictor
        vS = vS.oldTime()
           + dt*(1-offCentre)*dvSdt.oldTime() + dt*offCentre*dvSdt;
        //uS = TRiSK::ddirToFlux(vS);
        uS = Hops::ddirToFlux(vS);

        for (int corr=0; corr<nCorr; corr++)
        {
            #include "velocityVorticity.H"
            for(int icorr = 0; icorr < nNonOrthCorr; icorr++)
            {
                // explicit part of momentum change
                dvSdt = -pvFlux;
                if (!linear) dvSdt -= mesh.magSf()*fvc::snGrad(KE);
                
                vS = vS.oldTime()
                   + dt*(1-offCentre)*dvSdt.oldTime() + dt*offCentre*dvSdt;
                uS = Hops::ddirToFlux(vS)
                   - dt*offCentre*Hops::ddirToFluxOffDiag
                        (mesh.magSf()*fvc::snGrad(P));
                            
                fvScalarMatrix PEqn
                (
                    fvc::div(uS)
                  - fvm::laplacian(offCentre*dt*Hdiag, P)
                );
                PEqn.setReference(0, scalar(0));

                if (corr < nCorr-1) PEqn.solve();
                else PEqn.solve(mesh.solver(P.name() + "Final"));

                if (icorr == nNonOrthCorr-1) uS += PEqn.flux();
            }
            vS -= offCentre*dt*fvc::snGrad(P)*mesh.magSf();
        }

        #include "velocityVorticity.H"

        // update dvSdt for the next time-step
        dvSdt = -pvFlux -mesh.magSf()*fvc::snGrad(P);
        if (!linear) dvSdt -= mesh.magSf()*fvc::snGrad(KE);

        // velocity field for diagnostics
        Uf = Hops::faceReconstruct(uS);
        runTime.write();
//        #include "invariants.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
