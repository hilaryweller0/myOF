/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    addMountain

Description
    Read in a 3d mesh of the spherical atmosphere and add deform the mesh to
    create a mountain

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"

//using namespace Foam;
using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create spherical mesh for time = "
         << runTime.timeName() << Foam::nl << Foam::endl;

    fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

#   include "readEarthProperties.H"
    
    // Parameters for orography
    const dimensionedScalar earthRadius(earthProperties.lookup("earthRadius"));
    // Height at which mesh does not see orography
    const scalar meshOrogTop
         = readScalar(earthProperties.lookup("meshOrogTop"));
    const scalar mountainHeight
        = readScalar(earthProperties.lookup("mountainHeight"));
    const scalar mountainRadius
        = readScalar(earthProperties.lookup("mountainRadius_degrees"))*pi/180.;
    const scalar mountainLon
        = readScalar(earthProperties.lookup("mountainLon_degrees"))*pi/180.;
    const scalar mountainLat
        = readScalar(earthProperties.lookup("mountainLat_degrees"))*pi/180.;
    const vector omegaHat = unitVector(Omega.value());
    
    // New points and cell and face centres
    pointIOField newPoints
    (
        IOobject("points", mesh.time().timeName(), mesh.meshSubDir, mesh),
        mesh.points()
    );
    pointIOField newCellCentres
    (
        IOobject("cellCentres", mesh.time().timeName(), mesh.meshSubDir, mesh),
        mesh.C()
    );
    pointIOField newFaceCentres
    (
        IOobject("faceCentres", mesh.time().timeName(), mesh.meshSubDir, mesh),
        mesh.faceCentres()
    );
    
    // Move the points according to their distance from the ground and the 
    // orography
    
    forAll(newPoints, ip)
    {
        point& p = newPoints[ip];
        vector phat = unitVector(p);
        const scalar magp = mag(p);
        const scalar z = magp - earthRadius.value();
        const scalar plat = Foam::asin(vector(unitVector(p)) & omegaHat);
        const scalar plon = Foam::atan2
        (
            phat & vector(unitVector(vector(0.,1.,0.) - omegaHat.y()*omegaHat)),
            phat & vector(unitVector(vector(1.,0.,0.) - omegaHat.x()*omegaHat))
        );
        
        // Proportion of height wrt orography top
        const scalar zeta = 1 - z/meshOrogTop;
        
        // Orography height
        scalar h = 0;
        scalar dist = Foam::sqrt(sqr(plon - mountainLon) + sqr(plat - mountainLat));
        if (dist < mountainRadius)
        {
            h = mountainHeight*(1 - dist/mountainRadius);
        }
        
        if (zeta > 0)
        {
            // Proportion of orography applied at this altitude
            p *= (magp + zeta*h)/magp;
        }
    }
    
    newPoints.write();
    newCellCentres.write();
    newFaceCentres.write();
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
