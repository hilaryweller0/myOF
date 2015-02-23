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
    setT

Description
    Set T based on an array of temperatures for different heights. Linearly
    interpolate between the different heights

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    Info<< "Reading thermophysical properties\n" << endl;

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

    // The temperatures for different heights
    const scalarList TatHeight(envProperties.lookup("TatHeight"));
    // The heights for the temperatures
    const scalarList zT(envProperties.lookup("zT"));
    if (TatHeight.size() != zT.size())
    {
        FatalErrorIn("setT")
            << " size of TatHeight in environmentalProperties should be"
            << " the same as the size of zT"
            << exit(FatalError);
    }
        
    Info<< "Creating T\n" << endl;
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("T", dimTemperature, scalar(0))
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Loop over all cells and all boundaries and set T
    forAll(T, cellI)
    {
        const scalar z = -(mesh.C()[cellI] & ghat.value());
        // find this height within zT
        bool found = false;
        for (label il = 0; il < zT.size()-1 && !found; il++)
        {
            if (z >= zT[il] && z <= zT[il+1])
            {
                // linearly interpolate between T for levels il and il+1
                T[cellI] = TatHeight[il]*(zT[il+1]-z)/(zT[il+1]-zT[il])
                         + TatHeight[il+1]*(z-zT[il])/(zT[il+1]-zT[il]);
                found = true;
            }
        }
        if (!found)
        {
            FatalErrorIn("setT") << " cell at location "
                << mesh.C()[cellI] << " with height "
                << z << " not within range of given values " << zT
                << exit(FatalError);
        }
    }
    forAll(T.boundaryField(), patchI)
    {
        fvPatchField<scalar>& Tp = T.boundaryField()[patchI];
        forAll(Tp, facei)
        {
            const scalar z = -(mesh.C().boundaryField()[patchI][facei] & ghat.value());
            bool found = false;
            for (label il = 0; il < zT.size()-1 && !found; il++)
            {
                if (z >= zT[il] && z <= zT[il+1])
                {
                    // linearly interpolate between T for levels il and il+1
                    Tp[facei] = TatHeight[il]*(zT[il+1]-z)/(zT[il+1]-zT[il])
                             + TatHeight[il+1]*(z-zT[il])/(zT[il+1]-zT[il]);
                    found = true;
                }
            }
            if (!found)
            {
                FatalErrorIn("setT") << " boundary face at location "
                    << mesh.C().boundaryField()[patchI][facei]<<" with height "
                    << z << " not within range of given values " << zT
                    << exit(FatalError);
            }
        }
    }

    T.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

