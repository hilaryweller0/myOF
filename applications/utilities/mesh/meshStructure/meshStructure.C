/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Finds out if a mesh is structured (i.e. extruded) and writes maps if so.
    In parallel requires all columns of cells to be on a single processor.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "pointTopoDistanceData.H"
#include "PointEdgeWave.H"
#include "uindirectPrimitivePatch.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
tmp<GeoField> createField
(
    const labelIOList& lst,
    const typename GeoField::Mesh& mesh
)
{
    tmp<GeoField> tfld
    (
        new GeoField
        (
            IOobject
            (
                lst.name(),
                lst.time().timeName(),
                lst.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("a", dimless, 0)
        )
    );
    GeoField& fld = tfld();

    forAll(fld, i)
    {
        fld[i] = lst[i];
    }
    return tfld;
}


bool isStructuredCell
(
    const polyMesh& mesh,
    const labelList& faceToPatchEdgeAddressing,
    const labelList& pointLayer,
    const label layerI,
    const label cellI
)
{
    const cell& cFaces = mesh.cells()[cellI];

    // Count number of side faces
    label nSide = 0;
    forAll(cFaces, i)
    {
        if (faceToPatchEdgeAddressing[cFaces[i]] != -1)
        {
            nSide++;
        }
    }

    if (nSide != cFaces.size()-2)
    {
        return false;
    }

    // Check that side faces have correct point layers
    forAll(cFaces, i)
    {
        if (faceToPatchEdgeAddressing[cFaces[i]] != -1)
        {
            const face& f = mesh.faces()[cFaces[i]];

            label nLayer = 0;
            label nLayerPlus1 = 0;
            forAll(f, fp)
            {
                label pointI = f[fp];
                if (pointLayer[pointI] == layerI)
                {
                    nLayer++;
                }
                else if (pointLayer[pointI] == layerI+1)
                {
                    nLayerPlus1++;
                }
            }

            if (f.size() != 4 || (nLayer+nLayerPlus1 != 4))
            {
                return false;
            }
        }
    }

    return true;
}



// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "detect mesh layering"
    );

    #include "addRegionOption.H"
    argList::validArgs.append("(patches)");
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    #include "createNamedMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList patchIDs
    (
        pbm.patchSet(wordReList(IStringStream(args[1])())).sortedToc()
    );

    Info<< "Starting walk from patches "
        << UIndirectList<word>(pbm.names(), patchIDs)
        << nl
        << endl;

    label nFaces = 0;
    forAll(patchIDs, i)
    {
        nFaces += pbm[patchIDs[i]].size();
    }

    Info<< "Seeding " << returnReduce(nFaces, sumOp<label>()) << " patch faces"
        << nl << endl;


    labelList meshFaces(nFaces);
    nFaces = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = pbm[patchIDs[i]];
        forAll(pp, i)
        {
            meshFaces[nFaces++] = pp.start()+i;
        }
    }

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(mesh.faces(), meshFaces),
        mesh.points()
    );


    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Start of changes
    labelList patchFaces(pp.size());
    List<topoDistanceData> patchData(pp.size());
    forAll(pp, patchFaceI)
    {
        patchFaces[patchFaceI] = pp.addressing()[patchFaceI];
        patchData[patchFaceI] = topoDistanceData(patchFaceI, 0);
    }


    // Propagate information inwards
    FaceCellWave<topoDistanceData> distanceCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );


    // Write maps.

    labelIOList cellToPatchFaceAddressing
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nCells()
    );
    cellToPatchFaceAddressing.note() = "cell to patch face addressing";

    labelIOList cellLayer
    (
        IOobject
        (
            "cellLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nCells()
    );
    cellLayer.note() = "cell to layer addressing";

    forAll(cellToPatchFaceAddressing, cellI)
    {
        cellToPatchFaceAddressing[cellI] = cellData[cellI].data();
        cellLayer[cellI] = cellData[cellI].distance();
    }


    labelIOList faceToPatchFaceAddressing
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nFaces()
    );
    faceToPatchFaceAddressing.note() =
        "front/back face + turning index to patch face addressing";

    labelIOList faceToPatchEdgeAddressing
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nFaces()
    );
    faceToPatchEdgeAddressing.labelList::operator=(labelMin);
    faceToPatchEdgeAddressing.note() =
        "side face to patch edge addressing";

    labelIOList faceLayer
    (
        IOobject
        (
            "faceLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nFaces()
    );
    faceLayer.note() = "face to layer addressing";

    forAll(faceToPatchFaceAddressing, faceI)
    {
        label own = mesh.faceOwner()[faceI];
        if (mesh.isInternalFace(faceI))
        {
            label nei = mesh.faceNeighbour()[faceI];

            if (cellData[own].distance() == cellData[nei].distance())
            {
                // side face
                faceToPatchFaceAddressing[faceI] = 0;
                faceLayer[faceI] = cellData[own].distance();
            }
            else if (cellData[own].distance() < cellData[nei].distance())
            {
                // unturned face
                faceToPatchFaceAddressing[faceI] = faceI+1;
                faceToPatchEdgeAddressing[faceI] = -1;
                faceLayer[faceI] = faceData[faceI].distance();
            }
            else
            {
                // turned face
                faceToPatchFaceAddressing[faceI] = -(faceI+1);
                faceToPatchEdgeAddressing[faceI] = -1;
                faceLayer[faceI] = faceData[faceI].distance();
            }
        }
        else if (faceData[faceI].distance() == cellData[own].distance())
        {
            // starting face
            faceToPatchFaceAddressing[faceI] = -(faceI+1);
            faceToPatchEdgeAddressing[faceI] = -1;
            faceLayer[faceI] = faceData[faceI].distance();
        }
        else
        {
            // unturned face
            faceToPatchFaceAddressing[faceI] = faceI+1;
            faceToPatchEdgeAddressing[faceI] = -1;
            faceLayer[faceI] = faceData[faceI].distance();
        }
    }


    // Do all the faceToPatchEdgeAddressing that is still labelMin.
    labelIOList pointToPatchPointAddressing
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nPoints()
    );
    pointToPatchPointAddressing.note() =
        "point to patch point addressing";
    pointToPatchPointAddressing.labelList::operator=(-1);

    labelIOList pointLayer
    (
        IOobject
        (
            "pointLayer",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.nPoints()
    );
    pointLayer.note() = "point to layer addressing";



    {
        Info<< "Seeding " << returnReduce(pp.nPoints(), sumOp<label>())
            << " patch points"
            << nl << endl;

        // Field on edges and points.
        List<pointTopoDistanceData> edgeData(mesh.nEdges());
        List<pointTopoDistanceData> pointData(mesh.nPoints());

        // Start of changes
        labelList patchPoints(pp.nPoints());
        List<pointTopoDistanceData> patchData(pp.nPoints());
        forAll(pp.meshPoints(), patchPointI)
        {
            patchPoints[patchPointI] = pp.meshPoints()[patchPointI];
            patchData[patchPointI] = pointTopoDistanceData(patchPointI, 0);
        }


        // Walk
        PointEdgeWave<pointTopoDistanceData> distanceCalc
        (
            mesh,
            patchPoints,
            patchData,

            pointData,
            edgeData,
            mesh.globalData().nTotalPoints()  // max iterations
        );

        forAll(pointData, pointI)
        {
            pointToPatchPointAddressing[pointI] = pointData[pointI].data();
            pointLayer[pointI] = pointData[pointI].distance();
        }


        // Derive from originating patch points what the patch edges were.
        EdgeMap<label> pointsToEdge(pp.nEdges());
        forAll(pp.edges(), edgeI)
        {
            pointsToEdge.insert(pp.edges()[edgeI], edgeI);
        }

        // Look up on faces
        forAll(faceToPatchEdgeAddressing, faceI)
        {
            if (faceToPatchEdgeAddressing[faceI] == labelMin)
            {
                const face& f = mesh.faces()[faceI];

                // See if there is any edge
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    label nextPointI = f.nextLabel(fp);

                    EdgeMap<label>::const_iterator fnd = pointsToEdge.find
                    (
                        edge
                        (
                            pointData[pointI].data(), 
                            pointData[nextPointI].data()
                        )
                    );
                    if (fnd != pointsToEdge.end())
                    {
                        faceToPatchEdgeAddressing[faceI] = fnd();
                        // Note: could test whether the other edges on the
                        // face are consistent
                        break;
                    }
                }
            }
        }
    }


    // Use maps to find out mesh structure.
    {
        label nLayers = gMax(cellLayer)+1;
        labelListList layerToCells(invertOneToMany(nLayers, cellLayer));

        bool structured = true;
        forAll(layerToCells, layerI)
        {
            const labelList& lCells = layerToCells[layerI];

            forAll(lCells, lCellI)
            {
                label cellI = lCells[lCellI];

                structured = isStructuredCell
                (
                    mesh,
                    faceToPatchEdgeAddressing,
                    pointLayer,
                    layerI,
                    cellI
                );

                if (!structured)
                {
                    structured = false;
                    Info<< "This mesh is not structured starting at layer "
                        << layerI << nl << endl;
                    break;
                }
            }

            if (!structured)
            {
                break;
            }
        }

        if (structured)
        {
            Info<< "This mesh is structured with respect to patches "
                << UIndirectList<word>(pbm.names(), patchIDs) << nl << endl;
        }
    }



    Info<< "Writing maps for " << mesh.name()
        << " to " << mesh.facesInstance() << nl
        << endl;

    bool ok =
        cellToPatchFaceAddressing.write()
     && cellLayer.write()
     && faceToPatchFaceAddressing.write()
     && faceToPatchEdgeAddressing.write()
     && faceLayer.write()
     && pointToPatchPointAddressing.write()
     && pointLayer.write();

    if (!ok)
    {
        FatalErrorIn(args.executable())
            << "Failed writing maps for mesh " << mesh.name()
            << " at location " << mesh.facesInstance()
            << exit(FatalError);
    }



    //
    // Write maps as volFields,pointFields etc for postprocessing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //

    Info<< "Writing maps as volFields, surfaceFields, pointFields"
        << " for postprocessing" << nl
        << endl;

    {
        tmp<volScalarField> tcellMap
        (
            createField<volScalarField>(cellToPatchFaceAddressing, mesh)
        );
        tcellMap().correctBoundaryConditions();
        tcellMap().write();
    }
    createField<volScalarField>(cellLayer, mesh)().write();
    createField<surfaceScalarField>(faceToPatchFaceAddressing, mesh)().write();
    createField<surfaceScalarField>(faceToPatchEdgeAddressing, mesh)().write();
    createField<pointScalarField>
    (
        pointToPatchPointAddressing,
        pointMesh::New(mesh)
    )().write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
