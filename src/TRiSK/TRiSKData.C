/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "TRiSKData.H"
#include "fvc.H"
//#include "plane.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::TRiSKData::calcEllDist
(
    scalar& ell,
    scalar& dist,
    const label ie
) const
{
    const fvMeshWithDual& mesh = smesh();
    const polyPatch& bottomPatch = mesh.bottomPatch();

    // Points at either end of this edge
    const point p0 = unitVector
    (
        bottomPatch.localPoints()[bottomPatch.edges()[ie][0]]
    );
    const point p1 = unitVector
    (
        bottomPatch.localPoints()[bottomPatch.edges()[ie][1]]
    );
        
//    // Owner and neighbour cell centres
//    const point po = unitVector(bottomPatch.faceCentres()[edgeOwner(ie)]);
//    const point pn = unitVector(bottomPatch.faceCentres()[edgeNeib(ie)]);
    
    // Edge length
    ell = sphDist(p0, p1);

    // area associated with this edge
    scalar edgeArea = (mesh.edgeAreaOwn()[ie] + mesh.edgeAreaNei()[ie])
                     *0.25*sqr(mag(p0)+mag(p1));

    // d for this edge (distance between po and pn for conservation)
    dist = 2*edgeArea/ell;
    //dist = sphDist(po, pn);
}

void Foam::TRiSKData::calcStencilWeights
(
    labelList& stencil,
    scalarList& weight,
    const label faceI,
    const label cellI
) const
{
    const fvMeshWithDual& mesh = smesh();
//         = *(dynamic_cast<const fvMeshWithDual*>(&(this->mesh())));
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    const polyPatch& bottomPatch = mesh.bottomPatch();
    
    const scalarListList& facePointAreaFrac = mesh.facePointAreaFrac();
    //const scalarListList& facePointAreaFrac = facePointAreaFracs_;

    // elements on the patch that the cell and face refer to
    const label ie = mesh.faceToPatchEdge()[faceI];
    const label pf = mesh.cellToPatchFace()[cellI];
    //const label nSides = facePointAreaFracs_[pf].size();
    const label nSides = facePointAreaFrac[pf].size();

    // Circulate around patch face, starting from this edge and calc weights
    {
        const scalarList& areaFracsOld = facePointAreaFrac[pf];
        scalarList areaFracs(nSides);
        label ie0 = -1;
        for(label iee = 0; iee < nSides; iee++)
        {
            // first find this edge
            if (ie0 == -1 && bottomPatch.faceEdges()[pf][iee] == ie)
            {
                ie0 = iee;
            }
            // store the stencil and area fractions from this face onwards
            else if (ie0 >= 0)
            {
                stencil[iee-ie0-1] = bottomPatch.faceEdges()[pf][iee];
                areaFracs[iee-ie0] = areaFracsOld[iee];
            }
        }
        // and store the stencil and area fraction
        areaFracs[0] = areaFracsOld[ie0];
        for(label iee = 0; iee < ie0; iee++)
        {
            stencil[nSides-ie0-1+iee] = bottomPatch.faceEdges()[pf][iee];
            areaFracs[nSides-ie0+iee] = areaFracsOld[iee];
        }

        scalar areaFracSum = 0;
        for(label iee = 0; iee < nSides-1; iee++)
        {
            areaFracSum += areaFracs[iee];
            weight[iee] = areaFracSum - 0.5;
        }
    }
    
    // Map the patch stencil up to the mesh
    forAll(stencil, i)
    {
        bool faceFound = false;
        for(label fi = 0; fi < mesh.cells()[cellI].size() && !faceFound; fi++)
        {
            label faceI = mesh.cells()[cellI][fi];
            if (stencil[i] == mesh.faceToPatchEdge()[faceI])
            {
                stencil[i] = faceI;
                faceFound = true;
            }
        }
        if (!faceFound)
        {
            FatalErrorIn("TRiSKData::calcStencilWeights")
                << "cannot find face of cell " << cellI
                << " which maps onto one of the patch edges " << stencil
                << exit(FatalError);
        }
    }

    // Calculate the sign of the contribution from each face of the stencil
    
    // First calculate the direction indicator function for the owner faces
    const label tdirOwn = -sign
    (
        mesh.jdir()[faceI] & (mesh.rHatf()[stencil[0]] - mesh.rHatf()[faceI])
    )*sign(mesh.jdir()[faceI] & dualMesh.Sf()[mesh.dualFaceMap()[faceI]]);
    // Then calculate the sign of contribution for each owner face
    for(label iee = 0; iee < nSides-1; iee++)
    {
        label& ffsii = stencil[iee];
        label isOwner = mesh.faceOwner()[ffsii] == cellI ? 1 : -1;
        weight[iee] *= tdirOwn*isOwner;
    }
}



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TRiSKData, 0);
}


// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::TRiSKData::TRiSKData(const fvMesh& mesh__)
:
    MeshObject<fvMeshWithDual, Foam::MoveableMeshObject, TRiSKData>
        (*(dynamic_cast<const fvMeshWithDual*>(&mesh__))),
    smesh_(*(dynamic_cast<const fvMeshWithDual*>(&mesh__))),
    Htype_(NONE),
    ownerStencil_(mesh__.nFaces()),
    neighbourStencil_(mesh__.nInternalFaces()),
    vertToHorizStencil_(mesh__.nFaces(), FixedList<label, 4>(-1)),
    vertToHorizStencilSize_(mesh__.nFaces(), 0),
    own_weights_(mesh__.nFaces()),
    nei_weights_(mesh__.nInternalFaces()),
    circToFluxDiag_
    (
        IOobject("Hdiag", mesh__.time().timeName(), mesh__),
        mesh__,
        dimensionedScalar("H", dimless, scalar(1))
    ),
    circToFluxOffDiagWeights_(mesh__.nFaces()),
    circToFluxStencil_(mesh__.nFaces()),
    circToFluxStencilSize_(mesh__.nFaces()),
    ddirToFluxOffDiagWeights_(mesh__.nFaces()),
    ddirToFluxStencil_(mesh__.nFaces())
{
    const fvMeshWithDual& mesh = smesh();
    //const fvMeshWithDual& dualMesh = mesh.dualMesh();
    
    if (debug)
    {
        Info<< "Constructing TRiSKData from fvMesh"
            << endl;
    }

    // Calculate the owner and neighbour stencils for each vertical face
    forAll(ownerStencil_, faceI)
    {
        const label ie = mesh.faceToPatchEdge()[faceI];
        //const label faced = mesh.dualFaceMap()[faceI];
        if (ie != -1)
        {
            const label cellI = mesh.faceOwner()[faceI];
            const label cellSize = mesh.cells()[cellI].size();
            ownerStencil_[faceI].setSize(cellSize-3);
            own_weights_[faceI].setSize(cellSize-3);
            calcStencilWeights
            (
                ownerStencil_[faceI], own_weights_[faceI], faceI, cellI
            );
            
//            // Scale the weights by the area in the orthogonal direction
//            forAll(own_weights_[faceI], iee)
//            {
//                own_weights_[faceI][iee] /= mesh.Hdiag()[faceI];
//            }

//            // Scale the weights by l/d * dualSf/Sf
//            scalar dist = 0;
//            scalar ell = 0;
//            calcEllDist(ell, dist, ie);
//            const scalar distie = dist;
//            forAll(own_weights_[faceI], iee)
//            {
//                label eiee = mesh.faceToPatchEdge()[ownerStencil_[faceI][iee]];
//                calcEllDist(ell, dist, eiee);
//                own_weights_[faceI][iee] *= ell/distie
//                *dualMesh.magSf()[faced]/mesh.magSf()[ownerStencil_[faceI][iee]];
//            }
        }
    }
    forAll(neighbourStencil_, faceI)
    {
        const label ie = mesh.faceToPatchEdge()[faceI];
        //const label faced = mesh.dualFaceMap()[faceI];
        if (ie != -1)
        {
            const label cellI = mesh.faceNeighbour()[faceI];
            const label cellSize = mesh.cells()[cellI].size();
            neighbourStencil_[faceI].setSize(cellSize-3);
            nei_weights_[faceI].setSize(cellSize-3);
            calcStencilWeights
            (
                neighbourStencil_[faceI],nei_weights_[faceI],faceI, cellI
            );

//            // Scale the weights by the area in the orthogonal direction
//            forAll(own_weights_[faceI], iee)
//            {
//                nei_weights_[faceI][iee] /= mesh.Hdiag()[faceI];
//            }

//            // Scale the weights by l/d * dualSf/Sf
//            scalar dist = 0;
//            scalar ell = 0;
//            calcEllDist(ell, dist, ie);
//            const scalar distie = dist;
//            forAll(nei_weights_[faceI], iee)
//            {
//                label eiee = mesh.faceToPatchEdge()[neighbourStencil_[faceI][iee]];
//                calcEllDist(ell, dist, eiee);
//                nei_weights_[faceI][iee] *= ell/distie
//               *dualMesh.magSf()[faced]/mesh.magSf()[neighbourStencil_[faceI][iee]];
//            }
        }
    }
    
    // Assign the stencils for reconstructing w at vertical faces
    if (mesh.nGeometricD() == 3)
    {
        forAll(vertToHorizStencil_, faceI)
        {
            label ie = mesh.faceToPatchEdge()[faceI];
            if (ie != -1) // vertical face
            {
                label ivf = 0;
                const labelList& ownerCell = mesh.cells()[mesh.owner()[faceI]];
                forAll(ownerCell, iff)
                {
                    label faceII = ownerCell[iff];
                    if (mesh.faceToPatchEdge()[faceII] == -1)
                    {
                        vertToHorizStencil_[faceI][ivf++] = faceII;
                        vertToHorizStencilSize_[faceI]++;
                    }
                }
                const labelList& neibCell = mesh.cells()[mesh.neighbour()[faceI]];
                forAll(neibCell, iff)
                {
                    label faceII = neibCell[iff];
                    if (mesh.faceToPatchEdge()[faceII] == -1)
                    {
                        vertToHorizStencil_[faceI][ivf++] = faceII;
                        vertToHorizStencilSize_[faceI]++;
                    }
                }
            }
        }
    }
    
    // Read solution dict to determine type of non-orthogonal H vector
    const word Htypeword
    (
        mesh.solutionDict().lookupOrDefault<word>("nonOrthogHtype", "diagonal")
    );
    Htype_ = Htypeword == "diagonal" ? DIAGONAL :
             Htypeword == "Dubos"    ? DUBOS :
             Htypeword == "DubosWeller"    ? DUBOSW :
             Htypeword == "asymmetric" ? ASYMMETRIC : 
          //   Htypeword == "asymmetricCompact" ? ASYMMETRICCOMPACT:
             NONE;
    if (Htype_ == NONE)
    {
        FatalErrorIn("TRiSKData::TRiSKData")
            << "nonOrthogHtype in fvSolution must be one of diagonal, Dubos, "
            << "DubosWeller, asymmetric or asymmetricCompact, not " << Htypeword
            << exit(FatalError);
    }
    
    calcCircToFluxWeights();
}


bool Foam::TRiSKData::movePoints()
{
    return false;
}

