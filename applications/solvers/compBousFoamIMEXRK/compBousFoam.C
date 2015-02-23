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
    compBousFoam

Description
    Durran's non-linear, compressible Boussinesq equations in advective form
    solved using IMEX Runge Kutta time-stepping (explicit/implicit GWs,
    hevi/SI and linear/non-linear)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "TRiSK3d.H"
#include "butcherTableau.H"
#include "RKfield.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #define dt runTime.deltaT()
    #define aii imBt.coeffs()[iRK][iRK]
    
    const bool hevi(mesh.solutionDict().lookupOrDefault<bool>("hevi", false));
    const bool linear
    (
        mesh.solutionDict().lookupOrDefault<bool>("linear",false)
    );
    const bool SIgwaves
    (
        mesh.solutionDict().lookupOrDefault<bool>("SIgwaves", true)
    );
    const bool asymRK
    (
        mesh.solutionDict().lookupOrDefault<bool>("asymmetricRK", true)
    );
    
    #include "createFields.H"
    #include "createRKFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"
        
        for (int iRK = 0; iRK < exBt.nSteps(); iRK++)
        {
            // bouyancy equation (CP staggering so solving for bkSf on faces)
            bkSf = bkSf.oldTime() + dt*exBt.RKstep(iRK, sdBdt);
            if (SIgwaves) bkSf += dt*imBt.RKstep(iRK, fdBdt);
            
            streamFunc = streamFunc0*Foam::sin
            (
                omega*(runTime - dt*(1 -imBt.subTimes()[iRK]))
            );
            curlSf = fvc::interpolate(TRiSK3d::curl(streamFunc*vector(0,1,0)))
                   & mesh.Sf();

            fdudt[iRK] = bkSf + curlSf;
            
            // Explicit part of momentum equation for un
            un = G[iRK]*
            (
                un.oldTime()
              + dt*exBt.RKstep(iRK, sdudt)
              + dt*imBt.RKstep(iRK, fdudt)
              + dt*aii*fdudt[iRK]
            );
            
            // Continuity equation for P
            fvScalarMatrix PEqn
            (
                fvm::ddt(P)
              - exBt.RKstep(iRK, sdPdt)
              - imBt.RKstep(iRK, fdPdt)
              //+ cs2*aii*fvc::div(impDir*un)
              //+ cs2*aii*fvc::div(un)
              - fvm::laplacian(cs2*sqr(aii)*dt*impDir*G[iRK],P)
            );
            if (asymRK) PEqn += cs2*aii*fvc::div(un);
            else        PEqn += cs2*aii*fvc::div(impDir*un);
            
            PEqn.solve();
            
            // back substitute for un
            if (mag(aii) > SMALL) un += 1./(cs2*aii)*PEqn.flux();
            //un -= dt*aii*impDir*fvc::snGrad(P)*mesh.magSf();
            fdudt[iRK] -= impDir*fvc::snGrad(P)*mesh.magSf();
            
            // back substitute for bkSf
            if (SIgwaves)
            {
                fdBdt[iRK] = -un*N2*kdir;
                bkSf += dt*aii*fdBdt[iRK];
                fdudt[iRK] += dt*aii*fdBdt[iRK];
            }

            // Prognostic variables
            bk = fvc::reconstruct(bkSf);
            b = bk.component(vector::Z);
            U = fvc::reconstruct(un);
            Uf = fvc::interpolate(U);
            Uf += (un - (Uf&mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

            // update time derivatives
            sdBdt[iRK] = -kSf*(Uftot & TRiSK3d::reconstructVec(fvc::snGrad(b)));
            if (!SIgwaves) sdBdt[iRK] -= un*N2*kdir;

            sdudt[iRK] = -(Uftot & fvc::interpolate(fvc::grad(Uf))) & mesh.Sf();
            if (hevi) sdudt[iRK] -= expDir*fvc::snGrad(P)*mesh.magSf();
            if (hyperDiff)
            {
                sdudt[iRK] -= K*fvc::interpolate
                (
                    fvc::laplacian(delta2, fvc::laplacian(delta2, U))
                ) & mesh.Sf();
                sdBdt[iRK] -= K*fvc::interpolate
                (
                    fvc::laplacian(delta2, fvc::laplacian(delta2, b))
                )*kSf;
            }
    
            if (asymRK)
            {
                sdPdt[iRK] = -(Utot&fvc::grad(P));
                fdPdt[iRK] = - cs2*fvc::div(un);
            }
            else
            {
                sdPdt[iRK] = -(Utot&fvc::grad(P)) - cs2*fvc::div(expDir*un);
                fdPdt[iRK] = - cs2*fvc::div(impDir*un);
            }
        }
        
        divuNu = fvc::div(un);

        if (max(mag(imBt.weights())) > SMALL || max(mag(exBt.weights())) >SMALL)
        {
            // Final RK steps
            bkSf = bkSf.oldTime() + dt*exBt.RKfinal(sdBdt);
            if (SIgwaves) bkSf += dt*imBt.RKfinal(fdBdt);
            
            un = un.oldTime()
                + dt*exBt.RKfinal(sdudt)
                + dt*imBt.RKfinal(fdudt);
                
            P = P.oldTime() + dt*exBt.RKfinal(sdPdt) + dt*imBt.RKfinal(fdPdt);
            
            // Prognostic variables
            bk = fvc::reconstruct(bkSf);
            b = bk.component(vector::Z);
            U = fvc::reconstruct(un);
            Uf = fvc::interpolate(U);
            Uf += (un - (Uf&mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());

            divu = fvc::div(un);
        }
        else divu = divuNu;

        Info << "Sum b = " << (fvc::domainIntegrate(b)).value()/Vtot
             <<" Sum P = " << (fvc::domainIntegrate(P)).value()/Vtot
             <<" Sum U = " << (fvc::domainIntegrate(U)).value()/Vtot
             << endl;

        runTime.write();
                
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
