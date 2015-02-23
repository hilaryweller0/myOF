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
    polyDualPatch

Description
    Creates the dual of the specified patch, removes triangles and extrudes it. 

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
#include "removeTriangles.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMeshWithDual.H"

    // Identify the patch to dualise and then extrude
    const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);
    
    // Create the dual points and faces
    pointField dualPoints(ePatch.faceCentres());
    //pointField dualPoints(mesh.cellCentres());
    faceList   dualFaces(ePatch.nPoints());
    const labelListList& pointFaces = ePatch.pointFaces();
    
    forAll(dualFaces, ip)
    {
        const labelList& f = pointFaces[ip];
        dualFaces[ip].setSize(f.size());
        forAll(f, i)
        {
            dualFaces[ip][i] = f[i];
        }
        
        // Change order if necessary
        labelList& df = dualFaces[ip];
        for(bool noSwaps = false; !noSwaps;)
        {
            const point& a = dualPoints[df[0]];
            noSwaps = true;
            for(label i = 1; i < df.size()-1; i++)
            {
                const point& b = dualPoints[df[i]];
                const point& c = dualPoints[df[i+1]];
                // Flip direction if necessary
                if ((((b-a)^(c-b)) & (a+b+c)) < 0) 
                {
                    noSwaps = false;
                    label dfi = df[i];
                    df[i] = df[i+1];
                    df[i+1] = dfi;
                }
            }
        }
    }
    
    // Remove the triangles
    removeTriangles(dualFaces, dualPoints);
    
    PrimitivePatch<face, List, pointField> dualPatch(dualFaces, dualPoints);
    
    Info << "Creating extruded dual mesh" << endl;
    
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
    //const label nLayers = readLabel(earthProperties.lookup("nLayers"));
    
    extrudeModels::linearRadial radialExtrude(earthProperties);
    extrudedMesh dualMesh
    (
        IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime
        ),
        dualPatch,
        radialExtrude
    );
    
    // remove the empty zeroth boundary
    dualMesh.removeBoundary();
    DynamicList<polyPatch*> dualPatches(2);
    if (radialExtrude.nLayers() == 1)
    {
        dualPatches.append
        (
            new emptyPolyPatch
            (
                "originalPatch", dualPatch.size(),
                dualMesh.nInternalFaces(), 0, dualMesh.boundaryMesh(),
                "empty"
            )
        );
        dualPatches.append
        (
            new emptyPolyPatch
            (
                "otherSide", dualPatch.size(), 
                dualMesh.nInternalFaces() + dualPatch.size(),
                1,
                dualMesh.boundaryMesh(),
                "empty"
            )
        );
    }
    else
    {
        dualPatches.append
        (
            new polyPatch
            (
                "originalPatch", dualPatch.size(),
                dualMesh.nInternalFaces(), 0, dualMesh.boundaryMesh(),
                "patch"
            )
        );
        dualPatches.append
        (
            new polyPatch
            (
                "otherSide", dualPatch.size(), 
                dualMesh.nInternalFaces() + dualPatch.size(),
                1,
                dualMesh.boundaryMesh(),
                "patch"
            )
        );
    }

    dualMesh.addPatches(dualPatches);
    
    dualMesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
