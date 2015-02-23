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
    VoronoiSphereMesh

Description
    Voronoi mesher for the surface of a sphere. 
    Reads in either a list of points for the initial cell centres or cell
    centres from an existing OpenFOAM case.
    If necessary, performs various re-meshing tasks as directed in
    VoronoiSphereMeshDict.
    Renumbers the mesh using CGAL's spatial sorting
    Write out the primal mesh
    Create and write out the dual
    Write cell and face centres for primal and dual

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "InitialPoints.H"
#include "VoronoiSphereMeshing.H"
//#include "IOmanip.H"
//#include "extrudedMesh.H"
//#include "linearRadial.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    // Open control dictionary
    IOdictionary controlDict
    (
        IOobject
        (
            args.executable() + "Dict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    Info << "Reading initial points and required resolution" << endl;
    autoPtr<InitialPoints> initialPoints(InitialPoints::New(controlDict));

    Info << "Creating the Voronoi mesh from " << initialPoints().points().size()
         << " points" << endl;
    VoronoiSphereMeshing cgalMesh(controlDict, initialPoints);
    
    Info << "Adapting mesh resolution or optimisation" << endl;
    cgalMesh.adaptMesh();
    
    cgalMesh.write(runTime, Foam::fvMesh::defaultRegion, 0);
    //cgalMesh.write(runTime, "dualMesh", 1);
//    #include "createDualMesh.H"
//    #include "writeDualPrimalMaps.H"
//    #include "writeCellFaceCentres.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
