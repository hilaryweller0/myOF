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
    diamondMesh

Description
    Read in a mesh and decompose into diamonds centred around each face centre
    with points consisting of the face points and the cell centres either side

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "meshTools.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
//#include "OFstream.H"
//#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMeshWithDual.H"

    const scalar cellCentreRadius = mag(mesh.C()[0]);
    Info << "cellCentreRadius = " << cellCentreRadius << endl;

    // Identify the patch to diamondise and then extrude
    const polyPatch& pPatch(mesh.boundaryMesh()["originalPatch"]);
    const PrimitivePatch<face, Foam::List, pointField> ePatch
    (
        pPatch.localFaces(), pPatch.localPoints()
    );

    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    // Statistic of the old and new meshes
    const label nOldVerts = ePatch.nPoints();
    const label nOldCells   = ePatch.size();
    const label nOldEdges  = ePatch.nEdges();
    const label nDVerts = nOldVerts + nOldCells;
    const label nDCells = nOldEdges;

    pointField diamondPoints(nDVerts);
    faceList diamondFaces(nDCells);
    
    // New cell centres
    pointField diamondCentres(nDCells);

    // Set the diamond vertices
    for(label ip = 0; ip < nOldVerts; ip++)
    {
        diamondPoints[ip] = ePatch.points()[ip];
    }
    for(label ip = nOldVerts; ip < nOldVerts + nOldCells; ip++)
    {
        const label faceI = ip - nOldVerts;
        diamondPoints[ip] = ePatch.faceCentres()[faceI];
    }

    // Set the diamond faces
    for(label id = 0; id < nDCells; id++)
    {
        diamondFaces[id].setSize(4);

        diamondFaces[id][0] = ePatch.edges()[id][0];
        diamondFaces[id][1] = nOldVerts + ePatch.edgeFaces()[id][1];
        diamondFaces[id][2] = ePatch.edges()[id][1];
        diamondFaces[id][3] = nOldVerts + ePatch.edgeFaces()[id][0];
        
        diamondCentres[id] = cellCentreRadius*unitVector
        (
            ePatch.points()[ePatch.edges()[id][0]]
          + ePatch.points()[ePatch.edges()[id][1]]
        );
        
        // Check orientation
        point p1 = diamondPoints[diamondFaces[id][1]];
        vector e1 = p1 - diamondPoints[diamondFaces[id][0]];
        vector e2 = diamondPoints[diamondFaces[id][2]] - p1;
        if (((e1 ^ e2) & p1) < 0)
        {
            diamondFaces[id][3] = ePatch.edges()[id][0];
            diamondFaces[id][2] = nOldVerts + ePatch.edgeFaces()[id][1];
            diamondFaces[id][1] = ePatch.edges()[id][1];
            diamondFaces[id][0] = nOldVerts + ePatch.edgeFaces()[id][0];
        }
    }
    
    PrimitivePatch<face, List, pointField> newPatch(diamondFaces, diamondPoints);
    
    Info << "Creating new extruded mesh" << endl;
    
    extrudeModels::linearRadial radialExtrude(earthProperties);
    extrudedMesh newMesh
    (
        IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        ),
        newPatch,
        radialExtrude
    );

    #include "replaceBoundaries.H"

    Info << "Writing new extruded mesh" << endl;
    newMesh.write();
    
    #include "faceCellCentres.H"
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
