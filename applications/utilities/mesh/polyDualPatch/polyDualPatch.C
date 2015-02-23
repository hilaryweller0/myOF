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
    Creates the dual of the specified patch and extrudes it. 

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
#include "findEdgeMap.H"
#include "optimiseDual.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "optimiseDual",
        "number of iterations",
        "move dual points to reduce skewness"
    );

    argList::noParallel();
    #include "addRegionOption.H"
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    // Get times list
    instantList Times = runTime.times();
    #include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
    #include "createMeshWithDual.H"

//    enum CfType
//    {
//        ERROR,     // none of the following (0)
//        CROSS,     // primal-dual cross-over points (1)
//        PRIMAL,    // primal face centres (2)
//        DUAL       // dual face centres (3)
//    };

//    CfType Cftype = !args.optionFound("Cf") ? CROSS :
//                     args.optionRead<string>("Cf") == "cross" ? CROSS :
//                     args.optionRead<string>("Cf") == "primal" ? PRIMAL :
//                     args.optionRead<string>("Cf") == "dual" ? DUAL : ERROR;

    const label optimiseDual = !args.optionFound("optimiseDual") ? 0 :
                                args.optionRead<label>("optimiseDual");

    #include "createDualPatch.H"
    
    // Optimise dual mesh point locations to get cross-over points mid-way
    if (optimiseDual)
    {
        optimiseDualPoints(dualPoints, dualFaces, ePatch, optimiseDual);
    }

    PrimitivePatch<face, List, pointField> dualPatch(dualFaces, dualPoints);
    
    #include "createDualMesh.H"
    
    dualMesh.write();
    
    #include "dualPrimFaceMap.H"
    
    #include "faceCellCentres.H"
    
    //dualMesh.checkMesh(true);

    #include "cellPointMap.H"
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
