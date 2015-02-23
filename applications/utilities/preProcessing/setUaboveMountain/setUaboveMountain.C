/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    makeHotBubble

Description
    Modify the initial (uniform) theta to include a warm bubble as defined by
    Bryan and Fritsch, MWR 2002

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    static scalar piby2 = 0.5*constant::mathematical::pi;

    const surfaceVectorField& Cf = mesh.Cf();

    Info << "\nReading initialProperties" << endl;

    IOdictionary initialProperties
    (
        IOobject
        (
            "initialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar U0(readScalar(initialProperties.lookup("U0")));
    const scalar height1(readScalar(initialProperties.lookup("height1")));
    const scalar height2(readScalar(initialProperties.lookup("height2")));
    const scalar heightDiff = height2 - height1;

    Info<< "Reading Uf\n" << endl;
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );

    forAll(Uf, facei)
    {
        const scalar z = Cf[facei].z();
        if (z < height1) Uf[facei] = vector::zero;
        else if(z > height2) Uf[facei] = vector(U0, 0, 0);
        else
        {
            Uf[facei] = vector
            (
                sqr(Foam::sin(piby2*(z - height1)/heightDiff)), 0, 0
            );
        }
    }
    
    forAll(Uf.boundaryField(), patchI)
    {
        fvsPatchField<vector>& Ufp = Uf.boundaryField()[patchI];
        forAll(Ufp, facei)
        {
            const scalar z = Cf.boundaryField()[patchI][facei].z();
            if (z < height1) Ufp[facei] = vector::zero;
            else if(z > height2) Ufp[facei] = vector(U0, 0, 0);
            else
            {
                Ufp[facei] = vector
                (
                    sqr(Foam::sin(piby2*(z - height1)/heightDiff)), 0, 0
                );
            }
        }
    }

    Info<< "Writing Uf\n" << endl;
    Uf.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
