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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Global
    CourantNo

Description
    Calculates and outputs the maximum Courant Number.

\*---------------------------------------------------------------------------*/

    surfaceScalarField CoNum = mesh.surfaceInterpolation::deltaCoeffs()
                             *mag(phi/(hf*mesh.magSf()))*runTime.deltaT();

    Info<< "\nMax Courant Number = " << max(CoNum.internalField());

    // Gravity wave Courant number
    surfaceScalarField waveCoNum = mesh.surfaceInterpolation::deltaCoeffs()
                                 *sqrt(magg*hf)*runTime.deltaT();

    Info << "  Max gravity wave Courant number = "
         << max(waveCoNum.internalField()) << endl;

//scalar CoNum = 0.0;
//scalar meanCoNum = 0.0;

//if (mesh.nInternalFaces())
//{
//    scalarField sumPhi
//    (
//        fvc::surfaceSum(mag(phi))().internalField()/h.internalField()
//    );

//    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

//    meanCoNum =
//        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
//}

//Info<< "Courant Number mean: " << meanCoNum
//    << " max: " << CoNum << endl;

// ************************************************************************* //
