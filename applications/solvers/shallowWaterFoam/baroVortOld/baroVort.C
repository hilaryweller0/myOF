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
    nonOrthogSW

Description
    Transient Solver for multi-layer SW equation with non-orthogonal TRiSK

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK3d.H"
//#include "OFstream.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    #include "readSolutionDict.H"
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
            duSdt = -TRiSK3d::circToFlux(pvFlux)
                  + TRiSK3d::circToFluxOffDiag
                    (
                     dualMesh.magSf()*mesh.signedDualMap(fvc::snGrad(Phi))
                    );
            // - fvc::snGrad(KE);

            if (incompressible)
            {
                fvScalarMatrix PhiEqn
                (
                    fvc::div(duSdt)
                  - fvm::laplacian(mesh.Hdiag(), Phi)
                );
                PhiEqn.setReference(0, scalar(0));

                if (corr < nCorr-1) PhiEqn.solve();
                else PhiEqn.solve(mesh.solver(Phi.name() + "Final"));

                duSdt += PhiEqn.flux();
            }
            uS = uS.oldTime()
               + 0.5*runTime.deltaT()*(duSdt + duSdt.oldTime());

            // update full velocity (for diagnostics and for CLUST)
            Uf = TRiSK3d::faceReconstruct(uS/mesh.magSf(), mesh.idir());
            v = v.oldTime() - runTime.deltaT()*
            (
                pvFlux/dualMesh.magSf() - mesh.signedDualMap(fvc::snGrad(Phi))
            );

            if (!linear)
            {
                // Calculation of PV flux
                pv = f + TRiSK3d::curl(v);
                pvf = fvc::interpolate(pv);
            }
            //pvFlux = pvf*TRiSK3d::perp(uS);
            pvFlux = pvf*dualMesh.magSf()
              *TRiSK3d::speedMap(uS/mesh.magSf(), mesh.idir(), dualMesh.idir());
            if (conserveEnergy)
            {
                pvFlux = 0.5*
                (
                    pvFlux
                 // + TRiSK3d::perp(dualMesh.dualMap(pvf)*uS)
                 + dualMesh.magSf()*TRiSK3d::speedMap(dualMesh.dualMap(pvf)*uS/mesh.magSf(), mesh.idir(), dualMesh.idir())
                );
            }
        }
        
        uS = mesh.magSf()*TRiSK3d::speedMap(v, dualMesh.jdir(), mesh.idir());

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
