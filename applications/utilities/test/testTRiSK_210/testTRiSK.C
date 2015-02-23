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

#include "sphericalGeometry.H"
#include "TRiSK.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSphericalMesh.H"
    const TRiSKData& triskData = TRiSKData::New(mesh);
    const fvMesh& dualMesh = triskData.dualMesh();
    
//    Info << "primal points = " << mesh.points() << endl
//         << "primal C = " << mesh.C() << endl
//         << "primal Cf = " << mesh.Cf() << endl
//         << "dual points = " << dualMesh.points() << endl
//         << "dual C = " << dualMesh.C() << endl
//         << "dual Cf = " << dualMesh.Cf() << endl;
    
    const dimensionedScalar a("a", dimLength, mesh.bounds().mag());

    surfaceScalarField xf("xf", mesh.Cf().component(vector::X)/a);
    surfaceScalarField yf("yf", mesh.Cf().component(vector::Y)/a);
    surfaceScalarField zf("zf", mesh.Cf().component(vector::Z)/a);
    
    volScalarField y("y", dualMesh.C().component(vector::Y)/a);
    volScalarField z("z", dualMesh.C().component(vector::Z)/a);
    volScalarField r("r", mag(dualMesh.C())/a);

    Info << "Creating vector field Uf on faces" << endl;

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedVector("Uf", dimless, vector::zero)
    );
    if (max(mag(Uf)).value() < SMALL)
    {
//        Uf.replace(0, pow(yf+2*xf,3));
//        Uf.replace(1, pow(2*yf-xf,3));
//        Uf.replace(2, xf+yf+3*zf);
//        Uf -= (Uf & mesh.Cf())*mesh.Cf()/magSqr(mesh.Cf());

        Uf.replace(0, pow(yf+xf,2));
        Uf.replace(1, pow(yf-xf,2));
        Uf -= (Uf & mesh.Cf())*mesh.Cf()/magSqr(mesh.Cf());
    
        // Project the velocity to be divergence free
        volScalarField Phi
        (
            IOobject("Phi", runTime.timeName(), mesh),
            mesh, dimensionedScalar("Phi", dimLength, scalar(0)),
            "zeroGradient"
        );
        fvScalarMatrix pEqn(fvm::laplacian(Phi) == fvc::div(Uf & mesh.Sf()));
        pEqn.setReference(0, Phi[0]);
        pEqn.solve();
        Uf -= pEqn.flux()*mesh.Sf()/sqr(mesh.magSf());

        Uf.write();
    }

    // 2d stuff

    // Map the velocity onto the dual mesh
    surfaceVectorField Ufd("Uf", TRiSK::primalToDualFaceMap(Uf));
    Ufd.write();
    
    //Info << "Calculating the normal component of the velocity" << endl;
    surfaceScalarField u("u", Uf & mesh.Sf()/mesh.magSf());
    //u.write();
    
    Info << "Reconstructing velocity from the flux" << endl;
    surfaceVectorField Urecon
    (
        "Urecon",
        TRiSK::reconstructVec(u)
    );
    Urecon.write();

    Info << "Calculating the flux on the dual" << endl;
    surfaceScalarField phid("phi", Ufd & dualMesh.Sf());
    phid.write();
    
    Info << "Reconstructing the flux on the dual" << endl;
    surfaceScalarField phidr
    (
        "phir",
        dualMesh.magSf()*TRiSK::primalToDualFluxMap(TRiSK::perp(u))
    );
    phidr.write();
    
    Info << "Dual flux error" << endl;
    surfaceScalarField phidError("phiError", phidr - phid);
    phidError.write();
    
    Info << "Calculating divergence on the primal" << endl;
    volScalarField divu("divu", fvc::div(u*mesh.magSf()));
    divu.write();
    
    Info << "Mapping divergence to dual" << endl;
    volScalarField divuMapped("divMapped", TRiSK::primalToDualCellMap(divu));
    divuMapped.write();
    
    Info << "Calculating divergence on the dual from reconstructed flux" << endl;
    volScalarField divud("divud", fvc::div(phidr));
    divud.write();
    
    Info << "Dual divergence error" << endl;
    volScalarField divuError("divuError", divud - divuMapped);
    divuError.write();
    
    Info << "Curl of the velocity" << endl;
    volScalarField curlu("curlu", TRiSK::curl(u));
    curlu.write();
   
    Info << "KE" << endl;
    volScalarField KE("KE", TRiSK::ke(u));
    KE.write();
}

// ************************************************************************* //

