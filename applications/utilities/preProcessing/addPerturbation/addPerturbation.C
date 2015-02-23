// The FOAM Project // setJet.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   setJet
   \\  /           | Family: preProcessing
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
APPLICATION
    setJet

DESCRIPTION
    Add a pertuurbation to the height field, h

AUTHOR
    Hilary Weller

-------------------------------------------------------------------------------
*/

#include "fvMeshWithDual.H"
#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "zeroGradientPointPatchFields.H"
#include "slipFvPatchFields.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createMeshWithDual.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar degToRad = constant::mathematical::pi/180.;

    Info<< "Reading field h\n" << endl;
    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading earthProperties\n" << endl;
    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    
    const dimensionedScalar hhat(earthProperties.lookup("hhat"));
    const scalar lat2 = readScalar(earthProperties.lookup("lat2"))*degToRad;
    const scalar alpha = readScalar(earthProperties.lookup("alpha"));
    const scalar beta = readScalar(earthProperties.lookup("beta"));

    const volVectorField& C = mesh.C();
    
    // loop through all the cells and add the perturbation to the height
    forAll(C, celli)
    {
        dimensionedScalar a("a", dimLength, mag(C[celli]));
        scalar lat = Foam::asin(C[celli].z()/a.value());
        scalar lon = max(min(C[celli].y()/(a.value()*Foam::cos(lat)), 1.), -1);
        if(C[celli].x() > 0) lon = Foam::asin(lon);
        else lon = constant::mathematical::pi - Foam::asin(lon);

        h[celli] += hhat.value() * Foam::cos(lat) * Foam::exp(-sqr(lon/alpha))
                  * Foam::exp(-sqr((lat2 - lat)/beta));
    }
    
    // write out h
    h.rename("hPerturbed");
    h.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
