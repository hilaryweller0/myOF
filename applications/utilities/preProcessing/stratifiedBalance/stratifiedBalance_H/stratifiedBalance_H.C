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
    stratifiedBalance_dpdx

Description
    Find discretely balanced theta and Exner profiles given a set of ranges of
    N and then isothermal above the final N (Brunt Viaisally frequency)

\*---------------------------------------------------------------------------*/

#include "Hops.H"
#include "fvCFD.H"
#include "ExnerTheta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    Hops H(mesh);
    const surfaceScalarField gd("gd", g & H.delta());
    #include "createFields.H"
      
    const dictionary& itsDict = mesh.solutionDict().subDict("initialisation");
    const int maxIters = itsDict.lookupOrDefault<int>("maxIters", 100);
   
    // The Brunt Vaiasalla freqencies for the different layers
    const scalarList Brunt(envProperties.lookup("BruntVaisallaFreq"));
    // The extents of the different layers (size one greater)
    const scalarList zN(envProperties.lookup("zN"));
    if (Brunt.size()+1 != zN.size())
    {
        FatalErrorIn("setTheta")
            << " size of BruntVaisallaFreq in environmentalProperties should be"
            << " one smaller than the size of zN"
            << exit(FatalError);
    }
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        
    bool converged = false;
    for(label iter = 0; iter < maxIters && !converged; iter++)
    {
        Info << "Iteration " << iter << endl;
        
        // Loop over the different layers and set theta
        // values at the bottom of the next layer up
        scalar theta0 = T0.value();
        label cellBottom = -1;
        scalar zBottom = GREAT;
        for(label il = 0; il < Brunt.size(); il++)
        {
            forAll(theta, cellI)
            {
                const scalar z = mesh.C()[cellI].z();
                if (z >= zN[il] && z < zN[il+1])
                {
                    theta[cellI] = theta0*Foam::exp
                    (
                        sqr(Brunt[il])/mag(g.value())*(z-zN[il])
                    );
                }
                
                if (z > zN[il+1] && z < zBottom)
                {
                    zBottom = z;
                    cellBottom = cellI;
                }
            }
            
            forAll(theta.boundaryField(), patchI)
            {
                fvPatchField<scalar>& thetap = theta.boundaryField()[patchI];
                forAll(thetap, facei)
                {
                    const scalar z = mesh.C().boundaryField()[patchI][facei].z();
                    if (z >= zN[il] && z < zN[il+1])
                    {
                        thetap[facei] = theta0*Foam::exp
                        (
                            sqr(Brunt[il])/mag(g.value())*(z-zN[il])
                        );
                    }
                }
            }
            
            theta0 *= Foam::exp(sqr(Brunt[il])/mag(g.value())*(zN[il+1]-zN[il]));
            
//            Info << "After layer " << il << " cellBottom = " << cellBottom
//                 << " zBottom = " << zBottom << " theta0 = " << theta0
//                 << " ExnerBottom = " << Exner[cellBottom] << nl;
        }
        // Set theta to be isothermal above the top layer
        scalar Ttop = theta0*Exner[cellBottom];
        forAll(theta, cellI)
        {
            const scalar z = mesh.C()[cellI].z();
            if (z >= zBottom)
            {
                theta[cellI] = Ttop/Exner[cellI];
            }
            
            forAll(theta.boundaryField(), patchI)
            {
                fvPatchField<scalar>& thetap = theta.boundaryField()[patchI];
                forAll(thetap, facei)
                {
                    const scalar z = mesh.C().boundaryField()[patchI][facei].z();
                    if (z >= zBottom)
                    {
                        thetap[facei] = Ttop/Exner.boundaryField()[patchI][facei];
                    }
                }
            }
        }

        thetaf = fvc::interpolate(theta);
        rho = pRef/(R*theta)*pow(Exner, (1-kappa)/kappa);
        rhof = fvc::interpolate(rho);
        U = H.ddirToFlux(gd*rhof)
          - H.ddirToFluxOffDiag(Cp*rhof*thetaf*H.magd()*fvc::snGrad(Exner));
        fvScalarMatrix ExnerEqn
        (
            fvc::div(U)
          - fvm::laplacian(H.Hdiag()*Cp*rhof*thetaf*H.magd()/mesh.magSf(), Exner)
        );
        converged = ExnerEqn.solve(mesh.solver(Exner.name())).nIterations() == 0;
	}
                  
    Exner.write();
    theta.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

