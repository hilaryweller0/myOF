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
    Creates the dual of the specified patch and extrudes it

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sphericalGeometry.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
//#include "reorderMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createSphericalMesh.H"

//    // Write primal mesh cell and face centres for post processing
//    volVectorField C("C", mesh.C());
//    C.write();
//    surfaceVectorField Cf("Cf", mesh.Cf());
//    Cf.write();

    // Identify the patch to dualise and then extrude
    const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);

    // Create the dual points and faces
    pointField dualPoints(ePatch.faceCentres());
    //pointField dualPoints(mesh.cellCentres());
    faceList   dualFaces(ePatch.nPoints());
    const labelListList& pointFaces = ePatch.pointFaces();
    
//    // overwrite dualPoints with VoronoiPoints if they exist
//    IFstream vPts("VoronoiPoints");
//    if (vPts)
//    {
//        Info << "Reading VoronoiPoints" << endl;
//        vPts >> dualPoints;
//    }
//    
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
    
    // Remove the sides boundary and set the others to empty
    dualMesh.removeBoundary();
    DynamicList<polyPatch*> dualPatches(2);
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
    dualMesh.addPatches(dualPatches);
    
    dualMesh.write();
    
    // find map of faces from dual to primal mesh
    IOList<label> dualPrimalFaceMap
    (
        IOobject
        (
            "dualPrimalFaceMap",
            runTime.findInstance(dualMesh.meshDir(), "points"),
            dualMesh.meshDir(),
            runTime
        ),
        labelList(dualMesh.nInternalFaces(), -1)
    );
    
    for(label faceI = 0; faceI < dualMesh.nInternalFaces(); faceI++)
    {
        const label i = dualMesh.faceOwner()[faceI];
        const label j = dualMesh.faceNeighbour()[faceI];

        const labelList& iEdges = ePatch.pointEdges()[i];
        const labelList& jEdges = ePatch.pointEdges()[j];
        
        // In the primal mesh, find the edge between points i and j
        for(label iei = 0; iei < iEdges.size() && dualPrimalFaceMap[faceI] == -1; iei++)
        {
            for(label jej = 0; jej < jEdges.size() &&dualPrimalFaceMap[faceI]==-1;jej++)
            {
                if (iEdges[iei] == jEdges[jej])
                {
                    dualPrimalFaceMap[faceI] = iEdges[iei];
                }
            }
        }

        if(dualPrimalFaceMap[faceI]==-1 || dualPrimalFaceMap[faceI] >= dualMesh.nInternalFaces())
        {
            FatalErrorIn("polyDualPatch")
                << "Cannot find face in between primal mesh points "
                << i << " and " << j <<  " for dual face " << faceI << nl
                << "iEdges = " << iEdges << nl
                << "jEdges = " << jEdges << nl
                << exit(FatalError);
        }
    }
    dualPrimalFaceMap.write();
    
    // find map of faces from primal to dual mesh
    IOList<label> primalDualFaceMap
    (
        IOobject
        (
            "primalDualFaceMap",
            runTime.findInstance(dualMesh.meshDir(), "points"),
            dualMesh.meshDir(),
            runTime
        ),
        labelList(dualMesh.nInternalFaces(), -1)
    );
    
    forAll(dualPrimalFaceMap, i)
    {
        primalDualFaceMap[dualPrimalFaceMap[i]] = i;
    }
    primalDualFaceMap.write();

    Info << "Writing the dual mesh cell centres\n" << endl;
    
//    const scalar r1 = radialExtrude.Rsurface();
//    const scalar r2 = radialExtrude.Router();
    const scalar r1 = mag
    (
        radialExtrude(dualPatch.points()[0], dualPatch.pointNormals()[0], 0)
    );
    const label kOut = radialExtrude.nLayers();
    const scalar r2 = mag
    (
        radialExtrude(dualPatch.points()[0], dualPatch.pointNormals()[0], kOut)
    );
    //const scalar Rmid = Foam::sqrt((sqr(r1) + r1*r2 + sqr(r2))/3.);
    const scalar Rmid = 0.5*(r1 + r2);
    Info << "r1 = " << r1 << " r2 = " << r2 << " Rmid = " << Rmid << endl;

    pointIOField newCellCentres
    (
        IOobject("cellCentres", runTime.timeName(), "polyMesh", dualMesh),
        ePatch.localPoints()/mag(ePatch.localPoints())*Rmid
    );
    newCellCentres.write();

    Info << "Writing the dual mesh face centres\n" << endl;
    
    pointIOField newFaceCentres
    (
        IOobject("faceCentres", runTime.timeName(), "polyMesh", dualMesh),
        dualMesh.faceCentres()
    );
    for(label faceI = 0; faceI < dualMesh.nInternalFaces(); faceI++)
    {
        newFaceCentres[faceI] = mesh.faceCentres()[dualPrimalFaceMap[faceI]];
    }
    for(label cellI = 0; cellI < ePatch.nPoints(); cellI++)
    {
        label faceI = dualMesh.nInternalFaces() + cellI;
        const vector rHat = ePatch.localPoints()[cellI]/mag(ePatch.localPoints()[cellI]);
        newFaceCentres[faceI] = r1*rHat;
        newFaceCentres[faceI+ePatch.nPoints()] = r2*rHat;
    }
    newFaceCentres.write();
    
    // Check distance between dual and primal faces
    scalarList faceDist(mesh.nInternalFaces());
    forAll(faceDist, faceI)
    {
        faceDist[faceI] = mag
        (
            mesh.faceCentres()[faceI]
          - newFaceCentres[dualPrimalFaceMap[faceI]]
        );
    }
    Info << "Maximum distance between primal and dual face centres = "
         << max(faceDist)<< endl;
    
    //dualMesh.checkMesh(true);

    // Map from primal cell centres to dual points
    IOList<label> cellToPointMap
    (
        IOobject
        (
            "cellToPointMap",
            runTime.findInstance(dualMesh.meshDir(), "points"),
            dualMesh.meshDir(),
            runTime
        ),
        labelList(mesh.nCells(),-1)
    );

    for(label ip = 0; ip < mesh.nCells(); ip++)
    {
        cellToPointMap[ip] = dualPatch.whichPoint(ip%dualPatch.nPoints());
    }
    cellToPointMap.write();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
