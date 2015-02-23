/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "meshToPointField.H"
#include "centredCPCCellToCellStencilObject.H"
#include "upwindCPCCellToFaceStencilObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPointField<approxType>::meshToPointField
(
    const pointField& pts, const fvMesh& mesh, const stencilType sType
)
:
    mesh_(mesh),
    stencils_(),
    weights_(),
    stencilType_(sType)
{
    setPoints(pts);
}


template<class approxType>
Foam::meshToPointField<approxType>::meshToPointField
(
    const fvMesh& mesh, const stencilType sType
)
:
    mesh_(mesh),
    stencils_(),
    weights_(),
    stencilType_(sType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class approxType>
Foam::meshToPointField<approxType>::~meshToPointField()
{}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class approxType>
void Foam::meshToPointField<approxType>::setPointsCellCells
(
    const pointField& pts
)
{
    stencils_.setSize(pts.size());
    weights_.setSize(pts.size());

    forAll(pts, ip)
    {
        const point& p = pts[ip];
        label celli = mesh_.findCell(p);
        if (celli == -1)
        {
            celli = mesh_.findNearestCell(p);
        }
        label stencilSize = mesh_.cellCells()[celli].size()+1;
        stencils_[ip].setSize(stencilSize);
        weights_[ip].setSize(stencilSize);
    
        // Set the stencil to be the central cell plus its neighbours
        stencils_[ip][0] = celli;
        for(label is = 1; is < stencilSize; is++)
        {
            stencils_[ip][is] = mesh_.cellCells()[celli][is-1];
        }
    
        // Create the weights (depending on the approxType)
        approxType thisApproxType;
        thisApproxType.calcWeights(weights_[ip], mesh_, stencils_[ip], p);
    }
    
//    label ip = 2531;
//    Info << "Point " << ip << " has stencil " << stencils_[ip]
//         << " and weights " << weights_[ip] << endl;
}


template<class approxType>
void Foam::meshToPointField<approxType>::setPointsCellPointCells
(
    const pointField& pts
)
{
    stencils_.setSize(pts.size());
    weights_.setSize(pts.size());
    
    const extendedCentredCellToCellStencil& meshStencil
        = centredCPCCellToCellStencilObject::New(mesh_);

    // Write out the stencil size and order of accuracy for every face
    surfaceScalarField stencilSizeFeild
    (
        IOobject("stencilSize", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );
    surfaceScalarField orderField
    (
        IOobject("order", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );

    forAll(pts, ip)
    {
        const point& p = pts[ip];
        label celli = mesh_.findCell(p);
        if (celli == -1)
        {
            celli = mesh_.findNearestCell(p);
        }
        label stencilSize = meshStencil.stencil()[celli].size();
        label newStencilSize = 0;
        stencils_[ip].setSize(stencilSize);
        // Add cells only (not faces) to the stencils
        for(label is = 0; is < stencilSize; is++)
        {
            if (meshStencil.stencil()[celli][is] < mesh_.nCells())
            {
                stencils_[ip][newStencilSize++]
                     = meshStencil.stencil()[celli][is];
            }
        }
        stencilSize = newStencilSize;
        stencils_[ip].setSize(stencilSize);
        stencilSizeFeild[ip] = stencilSize;
    
        // Create the weights (depending on the approxType)
        weights_[ip].setSize(stencilSize);
        approxType thisApproxType;
        orderField[ip]
             = thisApproxType.calcWeights(weights_[ip], mesh_, stencils_[ip], p);
    }
    stencilSizeFeild.write();
    orderField.write();
}


template<class approxType>
void Foam::meshToPointField<approxType>::setPointsUpwindCPC
(
    const pointField& pts
)
{
    stencils_.setSize(pts.size());
    weights_.setSize(pts.size());
    
    const extendedUpwindCellToFaceStencil& meshStencil
        = upwindCPCCellToFaceStencilObject::New(mesh_, false, scalar(0.1));
    
    // Write out the stencil size and order of accuracy for every face
    surfaceScalarField stencilSizeFeild
    (
        IOobject("stencilSize", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );
    surfaceScalarField orderField
    (
        IOobject("order", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("s", dimless, scalar(1))
    );

    forAll(pts, ip)
    {
        const point& p = pts[ip];
        label celli = mesh_.findCell(p);
        if (celli == -1)
        {
            celli = mesh_.findNearestCell(p);
        }
        label faci = ip;
        
        // Work out whether to use the owner stencil or the neighbour stencil
        if (mesh_.owner()[faci] == celli)
        {
            stencils_[ip] = meshStencil.ownStencil()[faci];
        }
        else if (mesh_.neighbour()[faci] == celli)
        {
            stencils_[ip] = meshStencil.neiStencil()[faci];
        }
        else
        {
            FatalErrorIn("meshToPointField<approxType>::setPointsUpwindCPC")
                << "departure point " << ip << " at " << p
                << " not in owner cell or neighbour cell. Courant number too big"
                << exit(FatalError);
        }

        // Remove boundary faces from the stencil
        label newStencilSize = 0;
        for(label is = 0; is < stencils_[ip].size(); is++)
        {
            if (stencils_[ip][is] < mesh_.nCells())
            {
                stencils_[ip][newStencilSize++] = stencils_[ip][is];
            }
        }
        stencils_[ip].setSize(newStencilSize);
        
        stencilSizeFeild[ip] = newStencilSize;
    
        // Create the weights (depending on the approxType)
        weights_[ip].setSize(newStencilSize);
        approxType thisApproxType;
        orderField[ip]
             = thisApproxType.calcWeights(weights_[ip], mesh_, stencils_[ip], p);
    }
    stencilSizeFeild.write();
    orderField.write();
}

// ************************************************************************* //
