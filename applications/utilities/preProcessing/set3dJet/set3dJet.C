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
    set3dJet

Description
    Sets the velocity field for the Jablonowski Williamson steady jet

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"
//#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    
    #include "readInitProps.H"
    using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        u0*sqr(Foam::sin(2*mesh.latf()))*mesh.lonHatf()
    );
//    surfaceVectorField U("U", Uf);
//    U.write();
    Uf = TRiSK::divReconstruct(Uf & mesh.Sf());
    Uf.write();
    
//    surfaceVectorField Ufr
//    (
//        "Ufr",
//        TRiSK::faceReconstruct(Uf & mesh.idir(), mesh.idir())
//    );
//    Ufr.write();
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
