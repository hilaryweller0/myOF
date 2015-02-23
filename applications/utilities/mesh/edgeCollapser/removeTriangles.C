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
    removeTriangles

Description
    Remove triangles from a quasi-2d mesh (more triangles are created). 

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "removeFaces.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Get optional regionName
    word regionName;
    word regionDir;
    if (args.optionReadIfPresent("region", regionName))
    {
        regionDir = regionName;
        Info<< "Create mesh " << regionName << " for time = "
            << runTime.timeName() << nl << endl;
    }
    else
    {
        regionName = fvMesh::defaultRegion;
        Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;
    }
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    
    // List of faces to remove
    PackedBoolList removeFace(mesh.nFaces, false);
    
    const edgeList& edges = mesh.edges();
    //const pointField& points = mesh.points();

    PackedBoolList collapseEdge(mesh.nEdges());
    Map<point> collapsePointToLocation(mesh.nPoints());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        
        // if edge is connected to a face with 3 sides, collapse it
        collapseEdge[edgeI] = false;
        for
        (
            label fi = 0;
            fi < mesh.edgeFaces()[edgeI].size() && !collapseEdge[edgeI];
            fi++
        )
        {
            const label faceI = mesh.edgeFaces()[edgeI][fi];
            const face& edgeFace = mesh.faces()[faceI];
            if (edgeFace.size() == 3)
            //if (mesh.faceEdges()[faceI].size() == 3)
            {
                collapseEdge[edgeI] = true;
                
                // Collapse points to the triangle centre
                point triCentre = mesh.faceCentres()[mesh.edgeFaces()[edgeI][fi]];
                
                collapsePointToLocation.set(e[1], triCentre);
            }
        }
    }

    List<pointEdgeCollapse> allPointInfo;
    const globalIndex globalPoints(mesh.nPoints());
    labelList pointPriority(mesh.nPoints(), 0);

    collapser.consistentCollapse
    (
        globalPoints,
        pointPriority,
        collapsePointToLocation,
        collapseEdge,
        allPointInfo,
        true
    );

    // Topo change container
    polyTopoChange meshMod(mesh);

    // Put all modifications into meshMod
    bool anyChange = collapser.setRefinement(allPointInfo, meshMod);

    if (anyChange)
    {
        Info << "Writing mesh with collapsed edges" << endl;

        pointField samePoints = mesh.points();
        
        // Construct new mesh from polyTopoChange.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
        mesh.movePoints(samePoints);
        
        mesh.write();
    }    

    return 0;
}


// ************************************************************************* //
