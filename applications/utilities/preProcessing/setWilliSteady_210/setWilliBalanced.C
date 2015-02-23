// The FOAM Project // setWilliSteady.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   setWilliSteady
   \\  /           | Family: preProcessing
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
APPLICATION
    setWilliSteady

DESCRIPTION
    Set the initial field of U and h for puddleFoam.
    reads in and sets values for the 2nd Williamson et al shallow water test
    case. Global steady state non-linear zonal geostrophic flow

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/

#include "sphericalGeometry.H"
#include "TRiSK.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    const Switch overRideCellGeometry
    (
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "overRideCellGeometry", Switch(false)
        )
    );

    Foam::Info
        << "Create spherical mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvSphericalMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        overRideCellGeometry
    );

    #include "readSolutionDict.H"
    #include "readEarthProperties.H"
    const Switch sample(earthProperties.lookup("sampleInitialConditions"));

    const TRiSKData& triskData = TRiSKData::New(mesh);
    const fvSphericalMesh& dualMesh = triskData.dualMesh();
    const dimensionedScalar Rsphere("Rsphere", dimLength, mag(mesh.C()[0]));

    Info<< "Reading field h0\n" << endl;
    volScalarField h0
    (
        IOobject("h0", "constant", mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("h0", dimLength, 0.0)
    );

    const dimensionedScalar h1(earthProperties.lookup("h1"));
    const dimensionedScalar rotationPeriod
    (
        24*60*60*dimensionedScalar(earthProperties.lookup("rotationPeriod"))
    );
    const dimensionedScalar u0
    (
        mag(rotationPeriod.value()) > SMALL ?
            2*M_PI*Rsphere/rotationPeriod
          : dimensionedScalar("u0", dimVelocity, 0.)
    );
    const dimensionedScalar u0Bya = u0/Rsphere;
    const scalar au0 = Rsphere.value()*u0.value();
    const dimensionedScalar h2 = linear ?
        dimensionedScalar("h2", Rsphere*mag(Omega)*u0/magg) :
        dimensionedScalar("h2", (Rsphere*mag(Omega)*u0 +0.5*sqr(u0))/magg);
      
    Info << "h2 = " << h2 << "\nRsphere = " << Rsphere << endl;

    Info << "Calculating spherical geometry fields\n" << endl;

    const volVectorField rHat = mesh.C()/mag(mesh.C());
    const volVectorField rHatd = dualMesh.C()/mag(dualMesh.C());
    
    const volScalarField f = 2*(Omega & rHatd);

    Info<< "Creating field h\n" << endl;
    volScalarField h
    (
        IOobject("h", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("h", dimLength, 0.0)
    );
    
    Info<< "Creating field un\n" << endl;
    surfaceScalarField un
    (
        IOobject("un", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("un", dimVelocity, scalar(0))
    );
    
    Info<< "Creating field Uf\n" << endl;
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("Uf", dimVelocity, vector::zero)
    );
    surfaceVectorField UfD
    (
        IOobject("UfD", runTime.timeName(), dualMesh),
        dualMesh,
        dimensionedVector("UfD", dimVelocity, vector::zero)
    );

    if (!sample)
    {
        Info << "Creating field streamFunc\n" << endl;
        volScalarField streamFunc
        (
            IOobject("streamFunc", runTime.timeName(), dualMesh),
            dualMesh,
            dimensionedScalar("streamFunc", dimensionSet(0,2,-1,0,0), scalar(0))
        );

        Info << "Setting stream function on the dual\n" << endl;
        forAll(streamFunc, ip)
        {
            const scalar sinLat = rHatd[ip].z();
            streamFunc[ip] = -au0*sinLat;
        }
        streamFunc.write();
        const surfaceScalarField gradSf("gradSf", fvc::snGrad(streamFunc));
        gradSf.write();
        const surfaceVectorField Sfd("Sf", dualMesh.Sf());
        Sfd.write();

        Info << "Setting un\n" << endl;
        un = TRiSK::dualToPrimalFluxMap(gradSf, mesh);
        Uf = TRiSK::reconstructVec(un);
    }
    else
    {
        Info << "Sampling Uf on the faces\n" << endl;
        const surfaceVectorField& Cf = mesh.Cf();
        forAll(Uf, faceI)
        {
            scalar cosLat = Foam::sqrt(1 - sqr(Cf[faceI].z()/mag(Cf[faceI])));
            scalar u = u0.value()*cosLat;
            vector lonHat = unitVector(vector(0,0,1) ^ Cf[faceI]);
            Uf[faceI] = u*lonHat;
        }
        un = Uf & mesh.Sf()/mesh.magSf();
    }

    Info << "Setting internal values of h\n" << endl;
    
    forAll(h, cellI)
    {
        const scalar sinLat = rHat[cellI].z();
        //const scalar cosLat = Foam::sqrt(1 - sqr(sinLat));
        h[cellI] = (h1 - h2*sqr(sinLat)).value() - h0[cellI];
    }
    
    if (!sample)
    {
        // Make h in discrete geostrophic balance with un and div(h un) = 0
        //scalar eqnResidual = 1;
        //const scalar convergenceCriterion = 1e-13;

        volScalarField p
        (
            IOobject("p", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("p", dimensionSet(0,3,-1,0,0), scalar(0)),
            "zeroGradient"
        );

        bool converged = false;
        bool convergedOld = false;
        for(label it = 0; it < 30 && !convergedOld; it++)
        {
            Info << "\nIteration " << it << endl;
        
            surfaceScalarField hf = fvc::interpolate(h);
            if (linear) hf == H;
            surfaceScalarField hS = hf*mesh.magSf();

            volScalarField pv = linear ? f/H :
                (TRiSK::curl(un)+ f)/TRiSK::primalToDualCellMap(h);
            Uf = TRiSK::reconstructVec(un);
            UfD = TRiSK::primalToDualFaceMap(Uf);
            surfaceScalarField pvf
            (
                "pvf",
                TRiSK::dualToPrimalFaceMap(fvc::interpolate(pv),mesh)
            );
            surfaceScalarField pvFlux = conserveEnergy ? 
                    TRiSK::pvFlux(hf*un, pvf) : hf*pvf*TRiSK::perp(un);
            
            surfaceScalarField phidt = -hS*pvFlux;
            if (!linear) phidt -= hS*fvc::snGrad(TRiSK::ke(un));

            // Solve h equation to make h in geostrophic balance with un
            fvScalarMatrix hEqn
            (
                fvc::div(phidt)
              - fvm::laplacian(magg*hf, h)
              - fvc::laplacian(magg*hf, h0)
            );
            hEqn.setReference(0, h[0]);
            
            convergedOld = converged;
            converged = hEqn.solve(mesh.solver(h.name() + "Final")).nIterations()
                          == 0;
            //eqnResidual = hEqn.solve(mesh.solver(h.name() + "Final")).initialResidual();

            // Set up pressure equation to make h un divergence free
            surfaceScalarField phi = un*hS;
            fvScalarMatrix pEqn(fvm::laplacian(p) == fvc::div(phi));
            pEqn.setReference(0,p[0]);
            pEqn.solve(mesh.solver(p.name() + "Final"));
            phi -= pEqn.flux();
            un = phi/hS;
        }
    }
            
    h.write();
    un.write();
    Uf.write();

    Info << "Creating ff, Coriolis on the primal face centres" << endl;
    surfaceScalarField ff
    (
        IOobject("ff", runTime.constant(), mesh),
        2*(Omega & mesh.Cf())/mag(mesh.Cf())
    );
    ff.write();
    
    Info << "Creating ff, Coriolis on the dual face centres" << endl;
    surfaceScalarField ffd
    (
        IOobject("ff", runTime.constant(), dualMesh),
        //2*(Omega & dualMesh.Cf())/mag(dualMesh.Cf())
        TRiSK::primalToDualFaceMap(ff)
    );
    ffd.write();
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
