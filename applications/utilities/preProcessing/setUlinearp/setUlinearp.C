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
    setUlinearp

Description
    Set Uf based on an array of winds for different pressures. Linearly
    interpolate between the different pressures

\*---------------------------------------------------------------------------*/

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
    #include "createFields.H"
    
    // The temperatures for different heights
    const scalarList UatP(envProperties.lookup("UatP"));
    // The pressures for different winds
    const scalarList pU(envProperties.lookup("pU"));
    if (UatP.size() != pU.size())
    {
        FatalErrorIn("setUlinearp")
            << " size of UatP in environmentalProperties should be"
            << " the same as the size of pU"
            << exit(FatalError);
    }
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Loop over all faces and all boundaries and set Uf
    forAll(Uf, faceI)
    {
        const scalar pi =  pf[faceI];
        // find this pressure within pU
        bool found = false;
        for (label il = 0; il < pU.size()-1 && !found; il++)
        {
            if (sign((pi - pU[il])*(pi - pU[il+1])) <= SMALL)
            {
                // linearly interpolate between levels il and il+1
                Uf[faceI].x() = UatP[il]*(pU[il+1]-pi)/(pU[il+1]-pU[il])
                              + UatP[il+1]*(pi-pU[il])/(pU[il+1]-pU[il]);
                found = true;
            }
        }
        if (!found)
        {
            FatalErrorIn("setUlinearp") << " face with pressure "
                << pi << " not within range of given values " << pU
                << exit(FatalError);
        }
    }
    forAll(Uf.boundaryField(), patchI)
    {
        fvsPatchField<vector>& Up = Uf.boundaryField()[patchI];
        forAll(Up, facei)
        {
            const scalar pi = pf.boundaryField()[patchI][facei];
            bool found = false;
            for (label il = 0; il < pU.size()-1 && !found; il++)
            {
                if (sign((pi - pU[il])*(pi - pU[il+1])) <= SMALL)
                {
                    // linearly interpolate between levels il and il+1
                    Up[facei].x() = UatP[il]*(pU[il+1]-pi)/(pU[il+1]-pU[il])
                                  + UatP[il+1]*(pi-pU[il])/(pU[il+1]-pU[il]);
                    found = true;
                }
            }
            if (!found)
            {
                FatalErrorIn("setUlinearp") << " boundary face with pressure "
                    << pi << " not within range of given values " << pU
                    << exit(FatalError);
            }
        }
    }

    Uf.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

