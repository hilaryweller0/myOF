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
    Set the initial field of Uf and h for the 2nd Williamson et al shallow water
    test case. Global steady state non-linear zonal geostrophic flow

-------------------------------------------------------------------------------
*/

#include "meshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
//#include "mathematicalConstants.H"

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
        << "Create mesh with dual for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual mesh
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
    
    Foam::Info
        << "Create dualMesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual dualMesh
    (
        Foam::IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        mesh,
        overRideCellGeometry
    );
    mesh.setDual(dualMesh);
    
    const Switch linear
    (
        mesh.solutionDict().lookupOrDefault<Switch>("linear",false)
    );

    Info<< "Reading field h0\n" << endl;
    volScalarField h0
    (
        IOobject("h0", "constant", mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("h0", dimLength, 0.0)
    );

    Info << "Reading properties for initial conditions" << endl;
    const dimensionedScalar h1(mesh.earthProperties().lookup("h1"));
    const dimensionedScalar rotationPeriod
    (
        24*60*60*dimensionedScalar(mesh.earthProperties().lookup("rotationPeriod"))
    );
    const dimensionedScalar u0
    (
        mag(rotationPeriod.value()) > SMALL ?
            2*M_PI*mesh.earthRadius()/rotationPeriod
          : dimensionedScalar("u0", dimVelocity, 0.)
    );
    const dimensionedScalar u0Bya = u0/mesh.earthRadius();
    const dimensionedScalar h2 = linear ?
        dimensionedScalar
        (
            "h2", mesh.earthRadius()*mag(mesh.Omega())*u0/mesh.magg()
        ) :
        dimensionedScalar
        (
            "h2",
            (mesh.earthRadius()*mag(mesh.Omega())*u0 +0.5*sqr(u0))/mesh.magg()
        );

    Info << "Creating Uf" << endl;
    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        u0*Foam::cos(mesh.latf())*mesh.lonHatf()
    );
    Uf.write();

    Info << "Creating h" << endl;
    volScalarField h
    (
        IOobject("h", runTime.timeName(), mesh),
        h1 - h2*sqr(Foam::sin(mesh.lat())) - h0
    );
    h.write();
    
    Info << "Creating ff, Coriolis on the primal face centres" << endl;
    surfaceScalarField ff
    (
        IOobject("ff", runTime.timeName(), mesh),
        2*(mesh.Omega() & mesh.Cf())/mag(mesh.Cf())
    );
    ff.write();
    
    Info << "Creating ff, Coriolis on the dual face centres" << endl;
    surfaceScalarField ffd
    (
        IOobject("ff", runTime.timeName(), dualMesh),
        mesh.dualMap(ff)
        //2*(mesh.Omega() & dualMesh.Cf())/mag(dualMesh.Cf())
    );
    ffd.write();
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
