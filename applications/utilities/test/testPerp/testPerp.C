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
    testPerp

Description
    Tests the accuracy of the TRiSK perp operator

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK.H"
#include "Hoperator.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    
    Info << "Reading Uf" << endl;

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    
    // u in the normal direction for the primal and dual
    surfaceScalarField u("u", Uf & mesh.idir());
    surfaceScalarField ud("u", mesh.dualMap(Uf) & dualMesh.idir());
    
    // Compare with uperp
    surfaceScalarField uperp
    (
        "uperp",
        TRiSK::perp(u*mesh.magSf())/dualMesh.magSf()
    );

    Info << "Maximum perp error = " << max(mag(ud - uperp)).value() << endl;
    
    // Compare with uperp from reconstruct
    uperp = TRiSK::speedMap(u, mesh.idir(), dualMesh.idir());
    Info << "Maximum reconstruct perp error = " << max(mag(ud - uperp)).value()
         << endl;
    
//    // Recontruct Uf from normal components and write out
//    surfaceVectorField Ufr("Ufr", TRiSK::divReconstruct(u*mesh.magSf()));
//    Ufr.write();

    // Accuracy of H operator
    // v in the d direction (on the dual)
    surfaceScalarField vS
    (
        "vS", (mesh.dualMap(Uf) & dualMesh.jdir())*dualMesh.magSf()
    );
    // reconstruct u
    surfaceScalarField ur("ur", TRiSK::circToFlux(vS)/mesh.magSf());
    Info << "Maximum H error = " << max(mag(ur.internalField() - u.internalField()))
         << endl;
    Info << "Writing Herror " << endl;
    surfaceScalarField Herror("Herror", mag(ur - u));
    Herror.write();
    surfaceVectorField Ufr
    (
        "Ufr", ur*mesh.idir() + (Uf & mesh.jdir())*mesh.jdir()
    );
    Ufr.write();
}

// ************************************************************************* //

