// The FOAM Project // setResolution.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   setResolution
   \\  /           | Family: prrocessing
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
APPLICATION
    setJet

DESCRIPTION
    Set the required mesh density (resolution) as an analytic function for
    VoronoiSphereMesh

-------------------------------------------------------------------------------
*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading setResolutionDict" << endl;
    IOdictionary resDict
    (
        IOobject
        (
            "setResolutionDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // resolusion (dx in metres) in the fine and coarse regions
    const scalar minDx = readScalar(resDict.lookup("minDx"));
    const scalar maxDx = readScalar(resDict.lookup("maxDx"));
    
    // Centre of the fine region (in radians)
    const scalar lonCentre(readScalar(resDict.lookup("lonCentre"))*pi/180.);
    const scalar latCentre(readScalar(resDict.lookup("latCentre"))*pi/180.);
    
    // Inner and outer radius of the fine region (in metres)
    const scalar innerRadius(readScalar(resDict.lookup("innerRadius")));
    const scalar outerRadius(readScalar(resDict.lookup("outerRadius")));
    
    // the radius of the sphere
    const scalar a = mag(mesh.C()[0]);
    
    // the Cartesian form of the centre of the refined region
    const point rCentre
    (
        a*Foam::cos(lonCentre)*Foam::cos(latCentre),
        a*Foam::sin(lonCentre)*Foam::cos(latCentre),
        a*Foam::sin(latCentre)
    );
    
    Info<< "Creating requiredResolution" << endl;
    volScalarField rRes
    (
        IOobject
        (
            "requiredResolution",
            runTime.constant(),
            mesh
        ),
        mesh,
        dimensionedScalar("r", dimLength, maxDx)
    );

    // Loop through all the points and set fine resolution for the points in the
    // inner radius and a blend for those in the outer radius
    forAll(rRes, cellI)
    {
        // distance to the centre of the refined region
        scalar dist = mag(mesh.C()[cellI] - rCentre);
    
        if (dist < innerRadius)
        {
            rRes[cellI] = minDx;
        }
        else if (dist < outerRadius)
        {
            rRes[cellI] = ((dist-innerRadius)*maxDx + (outerRadius-dist)*minDx)
                    / (outerRadius - innerRadius);
        }
    }

    rRes.write();
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
