/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope tha it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    shallowWaterFoamVecInv

Description
    Transient solver for shallow water equations in vector invariant
    form using TRiSK operators or other operators

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sphericalGeometry.H"
#include "OFstream.H"
#include "TRiSK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSphericalMesh.H"
    #include "readSolutionDict.H"
    #include "readEarthProperties.H"
    const TRiSKData& triskData = TRiSKData::New(mesh);
    const fvMesh& dualMesh = triskData.dualMesh();
    #include "createFields.H"
    #include "initialInvariants.H"
    #include "invariants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    h == H*0;
    h[1] = 1;
    hv = TRiSK::primalToDualCellMap(h);
    Info //<< "phi = " << phi << endl
         //<< "KE = " << KE << endl
         //<< "f = " << f << endl
         << "hv = " << hv << endl
         << "pv = " << pv << endl
         << "pvFlux = " << pvFlux << endl;

    runTime++;
    hv.write();
    f.write();
    h.write();
    FatalErrorIn("") << exit(FatalError);

    Info<< "\nStarting time loop\n" << endl;
    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

#       include "CourantNo.H"

        // --- Outer loop

        for (int corr=0; corr<nCorr; corr++)
        {
            un = un.oldTime() - runTime.deltaT()*
            (
                offCentre*pvFlux + (1-offCentre)*pvFlux.oldTime()
              + (1-offCentre)*magg*fvc::snGrad(h.oldTime())
              + gsnGradh0
            );
            if (!linear)
            {
                un -= runTime.deltaT()*fvc::snGrad
                (
                    offCentre*KE + (1-offCentre)*KE.oldTime()
                );
            }

            phi = un*hS;
            
            // Helmholtz equations
            fvScalarMatrix hEqn
            (
                fvm::ddt(h)
              + fvc::div((1-offCentre)*phi.oldTime() + offCentre*phi)
              - fvm::laplacian(sqr(offCentre)*magg*hf*runTime.deltaT(), h)
            );

            if (corr < nCorr-1) hEqn.solve();
            else hEqn.solve(mesh.solver(h.name() + "Final"));

            // update phi and un
            phi += 1./offCentre*hEqn.flux();
            un = phi/hS;

            // update full velocity (for diagnostics and for CLUST)
            Uf = trisk ? TRiSK::reconstructVec(un)
                       : TRiSK::linReconstructVec(un);

            // Calculation of PV flux for trisk/not trisk, linear/non-linear
            // and conserveEnergy or not
            unD = trisk ? TRiSK::primalToDualFluxMap(TRiSK::perp(un))
                        : TRiSK::primalToDualFluxMap(TRiSK::linPerp(un));
            UfD = TRiSK::primalToDualFaceMap(Uf);

            if (!linear)
            {
                hv = TRiSK::primalToDualCellMap(h);
                KE = alphaKE*TRiSK::ke(un) + (1-alphaKE)*TRiSK::vertexKE(un);
//                hf = TRiSK::cellToFaceMap(h);
                hf = alphah*fvc::interpolate(h)
                   + (1-alphah)*TRiSK::dualToPrimalFaceMap
                     (fvc::interpolate(hv), mesh);
                hS = hf*mesh.magSf();
                
                if (!interpolateVorticity)
                {
                    pv = (TRiSK::curl(un) + f)/hv;
                    pvf = TRiSK::dualToPrimalFaceMap(fvc::interpolate(pv),mesh);
                }
                else
                {
                    vorticityByH = TRiSK::curl(un)/hv;
                    pv = vorticityByH + f/hv;
                    pvf = TRiSK::dualToPrimalFaceMap
                    (
                        fvc::interpolate(vorticityByH), mesh
                    ) + ff/hf;
                }
            }

            pvFlux = trisk ? pvf*TRiSK::perp(hf*un)
                           : pvf*TRiSK::linPerp(hf*un);
            if (conserveEnergy)
            {
                if (trisk) pvFlux = 0.5*(pvFlux + TRiSK::perp(pvf*hf*un));
                else       pvFlux = 0.5*(pvFlux + TRiSK::linPerp(pvf*hf*un));
            }
            
//            if (colocateU)
//            {
//                volVectorField U = fvc::reconstruct(un*mesh.magSf());
//                volVectorField pvp = unitVector(mesh.C())*TRiSK::primalCurl(un)
//                                   + 2*Omega;
//                volVectorField Coriolis = pvp ^ U;
//                Coriolis -= (Coriolis & unitVector(mesh.C()))*unitVector(mesh.C());
//                pvFlux = (fvc::interpolate(Coriolis) & mesh.Sf())/mesh.magSf();
//            }
        }

        #include "invariants.H"
        runTime.write();

        Info<< "\n    ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
            << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
