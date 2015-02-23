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
#include "TRiSK.H"
#include "OFstream.H"
#include "fvCFD.H"

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
    //#define HDIAG mesh.Hdiag()
    #define HDIAG trisk.circToFluxDiag()
    #include "createFields.H"
    #include "initialInvariants.H"
    if (errorDiags)
    {
        #include "invariants.H"
    }

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
                dvSdt = -pvFlux + gGradh0;
                if (!linear) 
                    dvSdt += mesh.dualFluxMap(fvc::snGrad(KE)*mesh.magSf());
                if (diffusion)
                {
                    del2u = -mesh.signedDualMap
                            (
                                fvc::snGrad(fvc::div(TRiSK::circToFlux(vS)))
                            );
                          - fvc::snGrad(TRiSK::curl(vS));
                
                    dvSdt += mu*dualMesh.magSf()*del2u;
                }
            
                vS = vS.oldTime()
                  + dt*(1-offCentre)*dvSdt.oldTime() + dt*offCentre*dvSdt;
            
                phi = TRiSK::circToFlux(hfd*vS)
                    + dt*offCentre*magg*TRiSK::circToFluxOffDiag
                      (
                          hfd*mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf())
                      );
                
                fvScalarMatrix hEqn
                (
                    fvm::ddt(h)
                  + fvc::div((1-offCentre)*phi.oldTime() + offCentre*phi)
                  - fvm::laplacian
                    (
                        sqr(offCentre)*magg*hf*dt*HDIAG,
                        h
                    )
                );

                if (corr < nCorr-1) hEqn.solve();
                else hEqn.solve(mesh.solver(h.name() + "Final"));

                if (icorr == nNonOrthCorr-1) phi += 1./offCentre*hEqn.flux();
            }

            vS += offCentre*dt*magg
                 *mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf());

            // update full velocity (for diagnostics and for CLUST)
            Ufd = TRiSK::faceReconstruct
            (
                vS/dualMesh.magSf(), dualMesh.jdir()
            );
            Ufd -= (Ufd & dualMesh.rHatf())*dualMesh.rHatf();
            Uf = dualMesh.dualMap(Ufd);
            //Uf = TRiSK::divReconstruct(phi/hf);
            
            if (!linear)
            {
                // Calculation of PV flux and KE
                hv = TRiSK::primalToDualCellMap(h);
                //hf = TRiSK::conserveInterp(h);
                hf = fvc::interpolate(h);
                hfd = mesh.dualMap(hf);
                phi = TRiSK::circToFlux(hfd*vS);
                phid = Ufd & dualMesh.Sf();
                KE = TRiSK::KE(TRiSK::circToFlux(vS), vS);
                if (!interpolateVorticity)
                {
                    pv = (f + TRiSK::curl(vS))/hv;
                    pvf = fvc::interpolate(pv);
                }
                else
                {
                    vorticityByH = TRiSK::curl(vS)/hv;
                    pv = vorticityByH + f/hv;
                    pvf = fvc::interpolate(vorticityByH) + ff/hfd;
                }
            }
            
            if (useTriskPerp)
            {
                pvFlux = pvf*TRiSK::perp(phi);
                if (conserveEnergy)
                {
                    pvFlux = 0.5*
                    (
                        pvFlux
                      + TRiSK::perp(dualMesh.dualMap(pvf)*phi)
                    );
                }
            }
            else
            {
                pvFlux = pvf*TRiSK::fluxMap(phi, mesh.idir(), dualMesh.idir());
            }
        }
        dvSdt = -pvFlux + gGradh0
              + magg*mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf());
        if (!linear)
            dvSdt += mesh.dualFluxMap(fvc::snGrad(KE)*mesh.magSf());
//        dvSdt += magg*mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf());
        if (diffusion)
        {
            del2u = -mesh.signedDualMap
                    (
                        fvc::snGrad(fvc::div(TRiSK::circToFlux(vS)))
                    );
                  - fvc::snGrad(TRiSK::curl(vS));
        
            dvSdt += mu*dualMesh.magSf()*del2u;
            del2U = TRiSK::faceReconstruct(del2u, dualMesh.jdir());
        }

        runTime.write();
        if (errorDiags)
        {
            #include "invariants.H"
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
