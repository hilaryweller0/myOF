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
    writeCellFaceCentres

Description
    Either writes cell centres and face centres as vol/surfaceVectorFields
    (writeGeoField option) or reads in vol/surfaceVectorFields cellCentres and
    faceCentres and writes them as vectorFields in polyMesh
    (writeMeshData option)
    If a file surfaceOnes exists the surface fields are multiplied by it

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("writeGeoField|writeMeshData");
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
    fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    word direction = args.args()[1];
    
    surfaceScalarField surfaceOnes
    (
        IOobject
        (
            "surfaceOnes",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar("one", dimless, scalar(1))
    );

    if (direction != "writeGeoField" && direction != "writeMeshData")
    {
        FatalErrorIn("writeCellFaceCentres")
            << "argument should either be writeGeoField or writeMeshData, not "
            << direction << exit(FatalError);
    }
    else if(direction == "writeGeoField")
    {
        Info<<"Writing cellCentres volVectorField from the cell centres"<<endl;
        volVectorField cellCentres
        (
            IOobject
            (
                "cellCentres",
                runTime.findInstance(mesh.meshDir(), "points"),
                mesh
            ),
            mesh,
            dimensionedVector("c", dimLength, vector::zero)
        );
        cellCentres.internalField() = mesh.cellCentres();
        cellCentres.write();

        Info<<"Writing faceCentres surfaceVectorField from the face centres"
            <<endl;
        surfaceVectorField faceCentres
        (
            IOobject
            (
                "faceCentres",
                runTime.findInstance(mesh.meshDir(), "points"),
                mesh
            ),
            mesh.Cf()*surfaceOnes,
            cellCentres.boundaryField().types()
        );
//        faceCentres.internalField() = mesh.faceCentres();
        faceCentres.write();
    }
    else // direction == writeMeshData
    {
        Info<<"Writing cellCentres vectorField from the cellCentres volVectorField" << endl;
        volVectorField cellCentres
        (
           IOobject("cellCentres",runTime.timeName(),mesh,IOobject::MUST_READ),
           mesh
        );
        pointIOField meshCellCentres
        (
            IOobject
            (
                "cellCentres",
                runTime.findInstance(mesh.meshDir(), "points"),
                mesh.meshSubDir,
                mesh
            ),
            cellCentres.internalField()
        );
        meshCellCentres.write();
        
        Info << "Writing faceCentres vectorField from the faceCentres surfaceVectorField" << endl;
        surfaceVectorField faceCentres
        (
           IOobject("faceCentres",runTime.timeName(),mesh,IOobject::MUST_READ),
           mesh
        );
        pointIOField meshFaceCentres
        (
            IOobject
            (
                "faceCentres",
                runTime.findInstance(mesh.meshDir(), "points"),
                mesh.meshSubDir,
                mesh
            ),
            mesh.nFaces()
        );
        for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            meshFaceCentres[faceI] = surfaceOnes[faceI]*faceCentres[faceI];
        }
        for(label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces();faceI++)
        {
            meshFaceCentres[faceI] = fieldAccess(faceCentres, faceI)
                                    *fieldAccess(surfaceOnes, faceI);
        }
        meshFaceCentres.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
