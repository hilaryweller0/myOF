/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "VoronoiSphereMeshing.H"
//#include "transform.H"
//#include "IFstream.H"
//#include "uint.H"
//#include "ulong.H"
////#include <list>
////#include "CGAL/squared_distance_3.h"
//#include "plane.H"
//#include "VectorSpaceFunctions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VoronoiSphereMeshing::VoronoiSphereMeshing
(
    const dictionary& dict,
    const InitialPoints& ips
)
:
    HTriangulation(),
    initialPoints(ips),
    radius(readScalar(dict.lookup("radiusOfSphere"))),
    writeIntermediate(dict.lookupOrDefault<Switch>("writeIntermediate", false)),
    VoronoiCentres(dict.lookupOrDefault<Switch>("VoronoiCentres", true)),
    relaxationFactorStart(dict.lookupOrDefault<scalar>("relaxationFactorStart",1)),
    relaxationFactorEnd(dict.lookupOrDefault<scalar>("relaxationFactorEnd",1)),
    doubleGlobalResolution(dict.lookupOrDefault<Switch>("doubleResolution", false)),
    addRemovePoints(dict.lookupOrDefault<Switch>("addRemovePoints", false)),
    changeTopology(dict.lookupOrDefault<Switch>("changeTopology", true)),
    fixBadShapes(dict.lookupOrDefault<Switch>("fixBadShapes", false)),
    nFakeLloyds(dict.lookupOrDefault<label>("nFakeLloyds", 0)),
    nLloydIterations(dict.lookupOrDefault<label>("nLloydIterations", 0)),
    nLaplacianSmooths(dict.lookupOrDefault<label>("nLaplacianSmooths", 0)),
    nHRiterations(dict.lookupOrDefault<label>("nHRIterations", 0)),
    maxPittewayIterations(dict.lookupOrDefault<label>("maxPittewayIterations", 0)),
    nTomitaSprings(dict.lookupOrDefault<label>("nTomitaSprings", 0)),
    beta(dict.lookupOrDefault<scalar>("TomitaSpringsBeta", 1.2)),
    smoothAll(dict.lookupOrDefault<Switch>("smoothAll", false)),
    spatialSort(dict.lookupOrDefault<Switch>("spatialSort", false)),
    minEdgeLength(dict.lookupOrDefault<scalar>("minEdgeLength", scalar(0))),
    meshOk()
{
    insertPoints(initialPoints.points());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VoronoiSphereMeshing::~VoronoiSphereMeshing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VoronoiSphereMeshing::adaptMesh
(
    Time& runTime,
    const word& regionName
)
{
    if (doubleGlobalResolution)
    {
        Info << "Doubling the resolution" << endl;
        doubleResolution();
    }

    if (addRemovePoints)
    {
        Info << "Adding and removing points to conform to the resolution" << endl;
        label nChanges = 0;
        nChanges += removePoints();
        nChanges += addPoints();
    }
    
//    if (maxEnlargeAnglesIters > 0)
//    {
//        Info << "Enlarge large angles in triangulation"
//             << endl;
//        enlargeAngles(maxEnlargeAnglesIters);
//    }
    
    if (nTomitaSprings) Info << "Tomita Spring relaxations" << endl;
    for (int iter=1; iter<=nTomitaSprings; iter++)
    {
        Info<<"  Iteration " << iter << ": ";
        TomitaSprings(relaxationFactor(iter, nTomitaSprings), changeTopology);
        if (fixBadShapes)
        {
            if (addRemovePoints) removeBadShapes();
            else deformBadShapes();
        }
    }

    if (nLaplacianSmooths) Info << "Laplacian smooths" << endl;
    for (int iter=1; iter<=nLaplacianSmooths; iter++)
    {
        Info<<"  Laplacian smooths " << iter << ": ";
        LaplacianSmooth();
        if (fixBadShapes)
        {
            if (addRemovePoints) removeBadShapes();
            else deformBadShapes();
        }
    }
        
    if (nFakeLloyds) Info << "Fake Lloyd iterations" << endl;
    for (int iter=1; iter<=nFakeLloyds; iter++)
    {
        Info<<"  Fake Lloyd iteration " << iter << ": ";
        LloydIteration(false);
    }
        
    if (nLloydIterations) Info << "Lloyd iterations" << endl;
    for (int iter=1; iter<=nLloydIterations; iter++)
    {
        Info<<"  Lloyd iteration " << iter << ": ";
        LloydIteration();
        if (fixBadShapes)
        {
            if (addRemovePoints) removeBadShapes();
            else deformBadShapes();
        }
        
        if (writeIntermediate)
        {
            runTime++;
            write(runTime, regionName, 0);
        }
    }
        
    if (nHRiterations) Info << "HR iterations" << endl;
    for (int iter=1; iter<=nHRiterations; iter++)
    {
        Info<<"  HR iteration " << iter << ": ";
        HRify();
        if (fixBadShapes)
        {
            if (addRemovePoints) removeBadShapes();
            else deformBadShapes();
        }
    }
    
    if (maxPittewayIterations) Info << "Pitteway iterations" << endl;
    label nNonPitt = 1;
    for (int iter = 1; iter <= maxPittewayIterations && nNonPitt > 0; iter++)
    {
        nNonPitt = PittewayIteration
        (
            relaxationFactor(iter, maxPittewayIterations)
        );
        if (fixBadShapes)
        {
            if (addRemovePoints) removeBadShapes();
            else deformBadShapes();
        }
    }
    
    if (spatialSort)
    {
        Info << "Renumber all the points using a CGAL spatial sort" << endl;
        spatialSortRenumber();
    }
}


// ************************************************************************* //
