/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Sets the mountain height for shallowWaterFoam.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMeshWithDual.H"
#include "volFields.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Foam::Info
        << "Create spherical mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        fvMeshWithDual::SPHERICALDIST
    );
    
    Info << "Creating mountain height, h0" << endl;
    volScalarField h0
    (
        IOobject("h0", runTime.findInstance("polyMesh", "points"), mesh),
        mesh,
        dimensionedScalar("h0", dimLength, scalar(0))
    );

    Info << "Reading earthProperties" << endl;
    IOdictionary earthProps
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );
    const dimensionedScalar hs0(earthProps.lookup("hs0"));
    const scalar radius(readScalar(earthProps.lookup("radius_degrees"))*pi/180.);
    const scalar lonC(readScalar(earthProps.lookup("lonC_degrees"))*pi/180.);
    const scalar latC(readScalar(earthProps.lookup("latC_degrees"))*pi/180.);
    
    forAll(h0, cellI)
    {
        scalar dist = Foam::sqrt
        (
            sqr(mesh.lon()[cellI] - lonC)
          + sqr(mesh.lat()[cellI] - latC)
        );
        
        if (dist < radius)
        {
            h0[cellI] = hs0.value()*(1 - dist/radius);
        }
    }
    
    h0.write();
        
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
