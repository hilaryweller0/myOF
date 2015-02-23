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

Application
    testMappedPatch

Description
    Test mapped b.c. by mapping face centres (mesh.C().boundaryField()).

\*---------------------------------------------------------------------------*/


#include "mappedPatchBase.H"
#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshTools.H"
#include "Time.H"
#include "OFstream.H"
#include "volFields.H"
#include "mappedPolyPatch.H"
#include "mappedFixedValueFvPatchFields.H"
#include "OBJstream.H"
#include "faceSet.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    //#include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Pout<< "Loaded " << mesh.name() << " with" << nl
        << "    cells   : " << mesh.nCells() << nl
        << "    faces   : " << mesh.nFaces() << nl
        << "    points  : " << mesh.nPoints() << nl
        << endl;

    fvMesh dualMesh
    (
        IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    Pout<< "Loaded " << dualMesh.name() << " with" << nl
        << "    cells   : " << dualMesh.nCells() << nl
        << "    faces   : " << dualMesh.nFaces() << nl
        << "    points  : " << dualMesh.nPoints() << nl
        << endl;


    const polyBoundaryMesh& primalPbm = mesh.boundaryMesh();
    label primalID = primalPbm.findPatchID("originalPatch");
    const polyPatch& primalPatch = primalPbm[primalID];

    const polyBoundaryMesh& dualPbm = dualMesh.boundaryMesh();
    label dualID = dualPbm.findPatchID("originalPatch");
    const polyPatch& dualPatch = dualPbm[dualID];


    // Construct a mapped patch 
    mappedPatchBase primalMapper
    (
        primalPatch,
        dualMesh.name(),                    //sampleRegion
        mappedPatchBase::NEARESTPATCHPOINT, //sampleMode
        dualPatch.name(),                   //samplePatch
        vector::zero                        //offset
    );
    mappedPatchBase dualMapper
    (
        dualPatch,
        mesh.name(),                        //sampleRegion
        mappedPatchBase::NEARESTPATCHPOINT, //sampleMode
        primalPatch.name(),                 //samplePatch
        vector::zero                        //offset
    );

    const mapDistribute& primalMap = primalMapper.map();

    // Use map to get point values onto faces.
    {
        pointField dualPatchPoints
        (
            dualPatch.points(),
            dualPatch.meshPoints()
        );

        // Get dual point locations in primal face order
        primalMap.distribute(dualPatchPoints);

        {
            OBJstream str(runTime.path()/"primalToDual.obj");

            forAll(primalPatch, faceI)
            {
                const point& fc = primalPatch.faceCentres()[faceI];
                const point& dualP = dualPatchPoints[faceI];


                Pout<< "Face:" << faceI << endl;
                Pout<< "    fc:" << fc << endl;
                Pout<< "    dualPatchPoints:" << dualP << endl;
                str.write(linePointRef(fc, dualP));
            }
        }
    }



    // Read all maps
    // ~~~~~~~~~~~~~
    const labelIOList cellToPatchFace
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList cellLayer
    (
        IOobject
        (
            "cellLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList faceToPatchFace
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList faceToPatchEdge
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList faceLayer
    (
        IOobject
        (
            "faceLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList pointToPatchPoint
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList pointLayer
    (
        IOobject
        (
            "pointLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );

    // Dual mesh
    const labelIOList dualCellToPatchFace
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualCellLayer
    (
        IOobject
        (
            "cellLayer",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualFaceToPatchFace
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualFaceToPatchEdge
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualFaceLayer
    (
        IOobject
        (
            "faceLayer",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualPointToPatchPoint
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );
    const labelIOList dualPointLayer
    (
        IOobject
        (
            "pointLayer",
            dualMesh.facesInstance(),
            dualMesh.meshSubDir,
            dualMesh,
            IOobject::MUST_READ
        )
    );


    // Get global point numbering for dual patch
    globalIndex globalDP(dualPatch.nPoints());
    labelList dualPointIDs(dualPatch.nPoints());
    forAll(dualPatch.meshPoints(), pointI)
    {
        dualPointIDs[pointI] = globalDP.toGlobal(pointI);
    }
    // Get dual (patch) point numbering in primal (patch) face order
    primalMap.distribute(dualPointIDs);


    // Get dual (patch) point numbering in primal cell order
    const labelList cellToDualPoint
    (
        UIndirectList<label>(dualPointIDs, cellToPatchFace)
    );

    // Swap 
    labelList nbrCellToDualPoint;
    syncTools::swapBoundaryCellList(mesh, cellToDualPoint, nbrCellToDualPoint);


//XXXXXX
    // In order of patch.faceEdges() construct dualPoint-dualPoints

    // For all the side faces on level0.

    edgeListList dualPointPoints(primalPatch.size());
    const labelListList& faceEdges = primalPatch.faceEdges();
    forAll(faceEdges, faceI)
    {
        dualPointPoints[faceI].setSize(faceEdges[faceI].size(), edge(-1, -1));
    }

    faceSet sideFaces0(mesh, "sideFaces0", 100);
    forAll(faceLayer, faceI)
    {
        if (faceLayer[faceI] == 0 && faceToPatchEdge[faceI] != -1)
        {
            sideFaces0.insert(faceI);

            label patchEdgeI = faceToPatchEdge[faceI];
            const labelList& patchFaces = primalPatch.edgeFaces()[patchEdgeI];

            edge dualEdge;
            label own = mesh.faceOwner()[faceI];
            if (mesh.isInternalFace(faceI))
            {
                label nei = mesh.faceNeighbour()[faceI];
                dualEdge = edge(cellToDualPoint[own], cellToDualPoint[nei]);
            }
            else
            {
                dualEdge = edge
                (
                    cellToDualPoint[own],
                    nbrCellToDualPoint[faceI-mesh.nInternalFaces()]
                );
            }
            forAll(patchFaces, i)
            {
                label patchFaceI = patchFaces[i];
                const labelList& fEdges = faceEdges[patchFaceI];
                label index = findIndex(fEdges, patchEdgeI);
                dualPointPoints[patchFaceI][index] = dualEdge;

                Pout<< "on patch face:" << primalPatch.faceCentres()[patchFaceI]
                    << " have connections to dual points "
                    << dualPointPoints[patchFaceI] << endl;
            }
        }

    }
    Info<< "Writing " << sideFaces0.size() << " faces to "
        << sideFaces0.objectPath() << endl;
    sideFaces0.write();


    // Reshuffle dualPointPoints into dual patch point order
    primalMap.reverseDistribute(dualPatch.nPoints(), dualPointPoints);

    // On the dual patch match edge to edge labels and send back
    EdgeMap<label> edgeToEdgeLabel(2*dualPatch.nEdges());
    globalIndex globalDE(dualPatch.nEdges());
    forAll(dualPatch.edges(), edgeI)
    {
        const edge& e = dualPatch.edges()[edgeI];
        edgeToEdgeLabel.insert
        (
            edge
            (
                globalDP.toGlobal(e[0]),
                globalDP.toGlobal(e[1])
            ),
            globalDE.toGlobal(edgeI)
        );
    }

    // Match dualPointPoints to the corresponding edge
    labelListList dualPointEdges(dualPointPoints.size());
    forAll(dualPointPoints, pointI)
    {
        const edgeList& pPoints = dualPointPoints[pointI];
        labelList& pEdges = dualPointEdges[pointI];

        pEdges.setSize(pPoints.size());
        forAll(pPoints, i)
        {
            pEdges[i] = edgeToEdgeLabel[pPoints[i]];
        }
    }

    // And send back to primal faces
    primalMap.distribute(dualPointEdges);
    // Since dualPointPoints was in order of faceEdges same goes for
    // dualPointEdges so now we've got the correspondence

    labelList primalEdgeToDualEdge(primalPatch.nEdges());
    forAll(faceEdges, faceI)
    {
        const labelList& fEdges = faceEdges[faceI];
        const labelList& dualPEdges = dualPointEdges[faceI];
        forAll(fEdges, i)
        {
            primalEdgeToDualEdge[fEdges[i]] = dualPEdges[i];
        }
    }


    List<Map<label> > compactMap;
    mapDistribute primalEdgeToDualEdgeMap
    (
        globalDE,
        primalEdgeToDualEdge,
        compactMap
    );


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
