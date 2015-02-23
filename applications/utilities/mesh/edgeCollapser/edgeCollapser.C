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
    edgeCollapser

Description
    Collapses edges which are shorter than minEdgeLength.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "edgeCollapser.H"
#include "fvMeshWithDual.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("minEdgeLength");
    #include "addRegionOption.H"
    #include "setRootCase.H"
    const scalar minEdgeLength = readScalar(IStringStream(args.args()[1])());
    #include "createTime.H"
    #include "createMeshWithDual.H"
    
    // Remember the old face centres
    pointField oldFaceCentres = mesh.faceCentres();
    
    // Edge collapsing engine
    edgeCollapser collapser(mesh);

    const edgeList& edges = mesh.edges();
    //const pointField& points = mesh.points();

    PackedBoolList collapseEdge(mesh.nEdges());
    Map<point> collapsePointToLocation(mesh.nPoints());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        
        // if edge is too short, collapse it
        collapseEdge[edgeI] = false;

        const scalar dist = e.mag(mesh.points());
        if (dist < minEdgeLength)
        {
            collapseEdge[edgeI] = true;
            collapsePointToLocation.set(e[1], mesh.points()[e[0]]);
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
        allPointInfo
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

        // Update fields
        mesh.updateMesh(map);

        mesh.movePoints(map().preMotionPoints());
        
        mesh.write();
        
        // create and write new face centres
        pointIOField faceCentres
        (
            IOobject("faceCentres", runTime.timeName(), "polyMesh", mesh),
            mesh.nFaces()
        );
        forAll(faceCentres, faceI)
        {
            faceCentres[faceI] = oldFaceCentres[map().faceMap()[faceI]];
        }
        faceCentres.write();
    }    

    return 0;
}


// ************************************************************************* //
