/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    testTRiSK

Description
    Tests the implementation of the TRiSK mesh and operators

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK3d.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
//    #include "createDualMesh.H"
    
    const dimensionedScalar a("a", dimLength, mesh.bounds().mag());

    surfaceScalarField xf("xf", mesh.Cf().component(vector::X)/a);
    surfaceScalarField yf("yf", mesh.Cf().component(vector::Y)/a);
    surfaceScalarField zf("zf", mesh.Cf().component(vector::Z)/a);

    Info << "Creating vector field Uf on faces" << endl;

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh), //, IOobject::MUST_READ),
        mesh,
        dimensionedVector("Uf", dimless, vector(1,1,0))
    );
//    Uf.replace(0, sqr(yf+xf));
//    Uf.replace(1, sqr(yf-xf));
    Uf.write();
    
    Info << "Reconstructing velocity vector from flux" << endl;
    surfaceVectorField Ufr
    (
        "Ufr",
        TRiSK3d::faceReconstruct(Uf & mesh.idir(), mesh.idir())
    );
    Ufr.write();

    Info << "Reconstructing velocity vector from circulation" << endl;
    surfaceVectorField Ufc
    (
        "Ufc",
        TRiSK3d::faceReconstruct(Uf & mesh.jdir(), mesh.jdir())
    );
    Ufc.write();
}

// ************************************************************************* //

