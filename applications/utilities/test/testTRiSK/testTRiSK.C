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
    const TRiSKData tRiSK3dData = TRiSKData::New(mesh);
    
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

//    // 3d stuff
//    // Map flux from idir to ddir and reconstruct velocity
//    surfaceScalarField uS("uS", (Uf & mesh.Sf()));
//    surfaceScalarField vS("vS", (Uf & mesh.ddir())*mesh.magSf());
////    surfaceVectorField Ufrecon("Ufrecon", TRiSK::divReconstruct3d(uS));
////    Ufrecon.write();

//    // Analytic and numerical curl
//    surfaceVectorField curl("curl3d", TRiSK::curl3d(vS));
//    curl.write();
//    surfaceVectorField curlA
//    (
//        IOobject("curlA", runTime.timeName(), mesh),
//        mesh,
//        dimensionedVector("curlA", curl.dimensions(), vector::zero)
//    );
//    curlA.replace(2, -4*yf/a);
//    curlA.write();
//    
//    volScalarField curlz("curlz", TRiSK::curl(mesh.dualFluxMap(vS)));
//    curlz.write();
//    volScalarField curlzA("curlzA", curlz);
//    curlzA = -4*y/a*z/r;
//    curlzA.write();
     

    // 2d stuff

    // Map the velocity onto the dual mesh and back
    surfaceVectorField Ufd("Uf", mesh.dualMap(Uf));
    Ufd.write();
    
//    surfaceVectorField Ufdd("Ufdd", dualMesh.dualMap(Ufd));
//    Ufdd.write();
    
    //Info << "Calculating the normal component of the velocity" << endl;
    surfaceScalarField u("u", Uf & mesh.idir());
    //u.write();
    
    Info << "u in the idir from v in the jdir" << endl;
    surfaceScalarField v("v", Ufd & dualMesh.jdir());
//    surfaceScalarField uji
//    (
//        "uji",
//        TRiSK::speedMap(v, dualMesh.jdir(), mesh.idir())
//    );
    //uji.write();
    surfaceScalarField uji2
    (
        "uji2",
        //TRiSK::circToFlux(v*dualMesh.magSf())/mesh.magSf()
        Hops::circToFlux(v*dualMesh.magSf())/mesh.magSf()
    );
//    uji2.write();
    
    Info << "Reconstructing velocity from the flux" << endl;
    surfaceVectorField Urecon
    (
        "Urecon",
        TRiSK::divReconstruct(Uf & mesh.Sf())
    );
    Urecon -= (Urecon & mesh.Cf())*mesh.Cf()/magSqr(mesh.Cf());
    Urecon.write();
    
    Info << "faceReconst of the velocity" << endl;
    surfaceVectorField Ufrecon
    (
        "Ufrecon",
        TRiSK::faceReconstruct(u, mesh.idir())
        //TRiSK::conserveInterp(fvc::reconstruct(u*mesh.magSf()))
        //linearInterpolate(fvc::reconstruct(u*mesh.magSf()))
    );
    Ufrecon -= (Ufrecon & mesh.Cf())*mesh.Cf()/magSqr(mesh.Cf());
    Ufrecon.write();
    
    Info << "faceReconst of the velocity using v" << endl;
    surfaceVectorField Uvrecon
    (
        "Uvrecon",
        TRiSK::faceReconstruct(v, dualMesh.jdir())
    );
    Uvrecon -= (Uvrecon & dualMesh.Cf())*dualMesh.Cf()/magSqr(dualMesh.Cf());
    Uvrecon.write();
    
    Info << "Reconstructing the velocity from uji2" << endl;
    surfaceVectorField Uji2recon
    (
        "Uji2recon",
        mesh.idir()*uji2 + mesh.jdir()*(Uf & mesh.jdir())
//        TRiSK::divReconstruct(uji2*mesh.magSf())
    );
    Uji2recon -= (Uji2recon & mesh.Cf())*mesh.Cf()/magSqr(mesh.Cf());
    Uji2recon.write();
    
    // Code to find the difference between different flux maps
    surfaceScalarField vS("vS", (Ufd & dualMesh.jdir())* dualMesh.magSf());

    Info << "Calculating the flux on the dual" << endl;
    surfaceScalarField phid("phi", Ufd & dualMesh.Sf());
    phid.write();
    
    Info << "Reconstructing the flux on the dual" << endl;
//    surfaceScalarField phidr("phir",  dualMesh.magSf()*TRiSK::perp(u));
    surfaceScalarField phidr("phir",  TRiSK::perp(u*mesh.magSf()));
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
    volScalarField curlu("curlu", TRiSK::curl(v*dualMesh.magSf()));
    curlu.write();
    
    Info << "uperp on the primal mesh" << endl;
    surfaceScalarField uperp
    (
        "uperp",
//        dualMesh.dualFluxMap(dualMesh.magSf()*TRiSK::perp(u))/mesh.magSf()
        dualMesh.dualFluxMap(TRiSK::perp(u*mesh.magSf()))/mesh.magSf()
    );
    uperp.write();
    
    Info << "divu on the dual from the linear reconstructed fluxes" << endl;
    surfaceScalarField dualFlux = mesh.dualMap(Urecon) & dualMesh.Sf();
    volScalarField divuLinRecon("divuLinRecon", fvc::div(dualFlux));
    divuLinRecon.write();
   
    Info << "KE" << endl;
    volScalarField KE("KE", TRiSK::KE(TRiSK::circToFlux(vS), vS));
    KE.write();
    
    Info << "gradDivu" << endl;
    surfaceVectorField gradDivu
    (
        "gradDivu",
        TRiSK::faceReconstruct
        (
           -mesh.signedDualMap
            (
                fvc::snGrad(fvc::div(TRiSK::circToFlux(vS)))
            ),
            dualMesh.jdir()
        )
    );
    gradDivu.write();
    
    surfaceVectorField curlCurlu
    (
        "curlCurlu",
        TRiSK::faceReconstruct
        (
            fvc::snGrad(TRiSK::curl(vS)),
            dualMesh.jdir()
        )
    );
    curlCurlu.write();
    
//    Info << "checking primalToDualCellMap" << endl;
//    volScalarField one
//    (
//        IOobject("one", runTime.timeName(), mesh),
//        mesh,
//        dimensionedScalar("one", dimless, scalar(1))
//    );
//    volScalarField oned = TRiSK::primalToDualCellMap(one);
//    Info << "1-oned goes from " << (min(1.-oned)).value()
//         << " to " << (max(1.-oned)).value()
//         << " and the global sum is " << sum(oned.internalField()*dualMesh.V())
//         << "\nwhereas the global sum of one is "
//         << sum(one.internalField()*mesh.V()) << endl;
//    
//    Info << setprecision(12)
//         << "primal mesh volume = " << sum(mesh.V()).value() << nl
//         << "dual mesh volume   = " << sum(dualMesh.V()).value() << nl
//         << "face vol sum       = " << sum(mesh.faceVol()).value() << nl
//         << "dual face vol sum  = " << sum(dualMesh.faceVol()).value() << endl;
//    
//    Info << "Checking conserveInterp" << endl;
//    surfaceScalarField onef = TRiSK::conserveInterp(one);
//    Info << "1-onef goes from " << (min(1.-onef)).value()
//         << " to " << (max(1.-onef)).value() << endl;
//    
//    Info << "Checking faceToCellMap" << endl;
//    volScalarField onefc = TRiSK::faceToCellMap(onef);
//    Info << "1-onefc goes from " << (min(1.-onefc)).value()
//         << " to " << (max(1.-onefc)).value() << endl;    
}

// ************************************************************************* //

