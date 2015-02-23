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
    createDiffuseTheta

Description
    Create a diffusion coefficient for theta (in the upper atmosphere)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    Info << "\nReading environmentalProperties" << endl;

    IOdictionary envProperties
    (
        IOobject
        (
            "environmentalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    dimensionedVector g(envProperties.lookup("g"));
    dimensionedVector ghat = g/mag(g);

    Info << "Reading in diffusion coefficient coefficients\n" << endl;
    const scalar zB(readScalar(envProperties.lookup("diffusionBase")));
    const scalar zt(readScalar(envProperties.lookup("diffusionTop")));
    const scalar diffBar(readScalar(envProperties.lookup("diffusionMean")));
        
    Info<< "Creating diffuseTheta\n" << endl;
    surfaceScalarField diffuseTheta
    (
        IOobject("diffuseTheta", runTime.constant(), mesh),
        mesh,
        dimensionedScalar("diffuseTheta", dimensionSet(1,-1,-1,0,0), scalar(0))
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Loop over all faces and set diffuseTheta
    forAll(diffuseTheta, faceI)
    {
        // First check if face does not have a vertical normal
        if(mag(mesh.Sf()[faceI] ^ ghat.value()) > mesh.magSf()[faceI]*1e-6)
        {
            // height of face centre
            const scalar z = -(mesh.Cf()[faceI] & ghat.value());
                
            // set the diffusion coefficient if the height is above base
            if (z > zB)
            {
                diffuseTheta[faceI] = diffBar*sqr(Foam::sin(0.5*pi*(z-zB)/(zt-zB)));
            }
            else if (z > zt)
            {
                FatalErrorIn("createDiffuseTheta") << "face " << faceI
                    << " has height " << z
                    << " but the diffusion is defined to lie between " << zB
                    << " and " << zt << exit(FatalError);
            }
        }
    }

    diffuseTheta.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

