/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    writes data in a lon lat z data table for the specified time step

Description
    

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("fieldName");
#   include "addTimeOptions.H"
    Foam::argList::addOption
    (
        "region",
        "meshRegion",
        "specify a non-default region to plot"
    );
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Foam::Info
        << "Create spherical mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        fvMeshWithDual::SPHERICALDIST
    );
    const word fieldName = args.args()[1];

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << " reading/writing field "
            << fieldName << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED
        );

        // initialise output file
        fileName outFile = args.rootPath() / args.caseName() / runTime.timeName();
        if (meshRegion != fvMesh::defaultRegion) outFile = outFile / meshRegion;
        outFile = outFile / fieldName + ".lonLatz";
        Info << "Writing file " << outFile << endl;
        OFstream os(outFile);
        os << "#lon    lat     z    " << fieldName << endl;


        if (fieldHeader.headerOk() && fieldHeader.headerClassName() == "volScalarField")
        {
            volScalarField vf(fieldHeader, mesh);

            for(label celli = 0; celli < mesh.nCells(); celli++)
            {
                os << mesh.lon()[celli] << "  " << mesh.lat()[celli] << "  "
                   << mesh.height()[celli] << "  " << vf[celli] << endl;
            }
        }
        else if (fieldHeader.headerClassName() == "surfaceVectorField")
        {
            surfaceVectorField vf(fieldHeader, mesh);

            for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                os << mesh.lonf()[faceI] << "  " << mesh.latf()[faceI] << "  "
                   << mesh.heightf()[faceI] << "  "
                   << (vf[faceI] & mesh.lonHatf()[faceI]) << "  "
                   << (vf[faceI] & mesh.latHatf()[faceI]) << "  "
                   << (vf[faceI] & mesh.rHatf()[faceI])
                   << endl;
            }
        }
        else if (fieldHeader.headerClassName() == "volVectorField")
        {
            volVectorField vf(fieldHeader, mesh);

            for(label cellI = 0; cellI < mesh.nCells(); cellI++)
            {
                os << mesh.lon()[cellI] << "  " << mesh.lat()[cellI] << "  "
                   << mesh.height()[cellI] << "  "
                   << (vf[cellI] & mesh.lonHat()[cellI]) << "  "
                   << (vf[cellI] & mesh.latHat()[cellI]) << "  "
                   << (vf[cellI] & mesh.rHat()[cellI])
                   << endl;
            }
        }
        else
        {
            Info << "No " << fieldName << endl;
        }
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
