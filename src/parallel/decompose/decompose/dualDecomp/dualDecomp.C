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

\*---------------------------------------------------------------------------*/

#include "dualDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"
#include "mappedPatchBase.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dualDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        dualDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualDecomp::dualDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict),
    methodDict_(decompositionDict_.subDict(typeName + "Coeffs")),
    region_(methodDict_.lookup("region")),
    patch_(methodDict_.lookup("patch"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dualDecomp::parallelAware() const
{
    return false;
}


Foam::labelList Foam::dualDecomp::decompose
(
    const polyMesh& dualMesh,
    const pointField& cc,
    const scalarField& cWeights
)
{
    const Time& runTime = dualMesh.time();

    // Loading primary region
    autoPtr<polyMesh> primalMeshPtr;
    
    if (!runTime.foundObject<polyMesh>(region_))
    {
        primalMeshPtr.reset
        (
            new polyMesh
            (
                Foam::IOobject
                (
                    region_,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    const polyMesh& primalMesh
    (
        primalMeshPtr.valid()
      ? primalMeshPtr()
      : runTime.lookupObject<polyMesh>(region_)
    );

    Info<< "Loaded primal mesh :" << region_
        << "  cells:" << primalMesh.globalData().nTotalCells()
        << "  faces:" << primalMesh.globalData().nTotalFaces()
        << "  points:" << primalMesh.globalData().nTotalPoints()
        << "  patches:" << primalMesh.boundaryMesh().size()
        << "  bb:" << primalMesh.bounds() << nl << endl;


    Info<< "Looking up patch " << patch_ << " on mesh " << primalMesh.name()
        << " and " << dualMesh.name() << nl << endl;

    const polyPatch& primalPatch = primalMesh.boundaryMesh()[patch_];
    const polyPatch& dualPatch = dualMesh.boundaryMesh()[patch_];


    Info<< "Loading extrusion maps" << nl << endl;

    const labelIOList primalDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            primalMesh.facesInstance(),
            primalMesh,
            IOobject::MUST_READ
        )
    );

    const labelIOList dualPointToPatch
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


    // Get mapping from primal to dual
    Info<< "Constructing primal to dual patch mapper" << nl << endl;

    mappedPatchBase primalMapper
    (
        primalPatch,
        dualMesh.name(),                    //sampleRegion
        mappedPatchBase::NEARESTPATCHPOINT, //sampleMode
        patch_,                             //samplePatch
        vector::zero                        //offset
    );
    const mapDistribute& primalMap = primalMapper.map();


    // Get decomposition of patch faces
    labelList decomp(primalDecomposition, primalPatch.faceCells());

    // Push primal decomposition onto dual points
    primalMap.reverseDistribute(dualPatch.nPoints(), decomp);


    // Decompose dual points accordingly
    labelList finalDecomp(dualMesh.nCells());
    forAll(dualMesh.cells(), cellI)
    {
        const cell& cFaces = dualMesh.cells()[cellI];

        label minProcI = labelMax;

        forAll(cFaces, cFaceI)
        {
            const face& f = dualMesh.faces()[cFaces[cFaceI]];

            forAll(f, fp)
            {
                label patchPointI = dualPointToPatch[f[fp]];
                minProcI = min(minProcI, decomp[patchPointI]);
            }
        }
        finalDecomp[cellI] = minProcI;
    }
    
    return finalDecomp;
}


Foam::labelList Foam::dualDecomp::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pointWeights
)
{
    notImplemented
    (
        "dualDecomp::decompose\n"
        "(\n"
        "    const labelListList&,\n"
        "    const pointField&,\n"
        "    const scalarField&\n"
        ")\n"
    );

    return labelList::null();
}


// ************************************************************************* //
