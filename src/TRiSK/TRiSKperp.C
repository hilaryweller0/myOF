/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "TRiSKperp.H"
#include "TRiSKData.H"
#include "TRiSKFaceToCellMap.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceScalarField> perp(const surfaceScalarField& uS)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    const TRiSKData& TRiSKData = TRiSKData::New(mesh);
    const fvMeshWithDual& dualMesh = mesh.dualMesh();

    tmp<surfaceScalarField> treconField
    (
        new surfaceScalarField
        (
            IOobject(uS.name(), uS.instance(), dualMesh),
            dualMesh,
            dimensionedScalar(uS.name(), uS.dimensions(), scalar(0)),
            uS.boundaryField().types()
        )
    );
    surfaceScalarField& uSd = treconField();
    
    forAll(uSd, faced)
    {
        const label faceI = dualMesh.dualFaceMap()[faced];
        if (faceI != -1)
        {
            // Contributions from primal owner cell faces
            const scalarList& weights = TRiSKData.own_weights()[faceI];
            const labelList& stencil = TRiSKData.ownerStencil()[faceI];
        
            forAll(stencil, i)
            {
                uSd[faced] += fieldAccess(uS, stencil[i])*weights[i];
            }
        }
    }
    forAll(uSd.boundaryField(), patchI)
    {
        label faced = dualMesh.boundaryMesh()[patchI].start()-1;
        
        forAll(uSd.boundaryField()[patchI], patchFace)
        {
            faced++;
            const label faceI = dualMesh.dualFaceMap()[faced];

            if (faceI != -1)
            {
                const scalarList& weights = TRiSKData.own_weights()[faceI];
                const labelList& stencil = TRiSKData.ownerStencil()[faceI];
        
                forAll(stencil, i)
                {
                    uSd.boundaryField()[patchI][patchFace]
                        += fieldAccess(uS, stencil[i])*weights[i];
                }
            }
        }
    }

    // Contributions from primal neighbour cell faces
    forAll(uSd, faced)
    {
        const label faceI = dualMesh.dualFaceMap()[faced];
        
        if (faceI != -1)
        {
            const scalarList& weights = TRiSKData.nei_weights()[faceI];
            const labelList& stencil  = TRiSKData.neighbourStencil()[faceI];
            
            forAll(stencil, i)
            {
                uSd[faced] += fieldAccess(uS, stencil[i])*weights[i];
            }
        }
    }

    return treconField;
}


tmp<surfaceScalarField> perp(const tmp<surfaceScalarField>& tuS)
{
    tmp<surfaceScalarField> tuSperp
    (
        TRiSK::perp(tuS())
    );
    tuS.clear();
    return tuSperp;
}


tmp<surfaceVectorField> perpVec(const surfaceScalarField& uS)
{
    FatalErrorIn("TRiSK::perpVec")
        << "This version has errors. See divReconstruct3d for correct version"
        << abort(FatalError);

    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));

    const TRiSKData& TRiSKData = TRiSKData::New(mesh);
    const fvMeshWithDual& dualMesh = mesh.dualMesh();

    // horizontal contribution to vertical faces
    tmp<surfaceVectorField> treconField
    (
        new surfaceVectorField
        (
            IOobject(uS.name(), uS.instance(), mesh),
            uS*mesh.idir()/mesh.magSf()
          + dualMesh.dualMap(perp(uS)*dualMesh.idir()/dualMesh.magSf())
        )
    );
    surfaceVectorField& Uf = treconField();
    
    // vertical contribution to vertical faces
    forAll(Uf, faceI)
    {
        //const label faceI = dualMesh.dualFaceMap()[faced];
        //if (faceI != -1)
        const label stencilSize =TRiSKData.vertToHorizStencilSize()[faceI];
        if
        (
            mesh.faceToPatchEdge()[faceI] != -1
         && stencilSize != 0
        )
        {
            const FixedList<label, 4>& stencil
                 = TRiSKData.vertToHorizStencil()[faceI];
            scalar w = 0;
            for(label is = 0; is < stencilSize; is++)
            {
                w += fieldAccess(uS, stencil[is]);
            }
            w /= scalar(stencilSize);
        
            Uf[faceI] += w*mesh.kdir()[faceI];
        }
    }

    volVectorField U = faceToCellMap(Uf);
    
    // Uf for horizontal faces
    forAll(Uf, faceI)
    {
        if (mesh.faceToPatchEdge()[faceI] == -1)
        {
            Uf[faceI] += 0.5*
                (U[mesh.owner()[faceI]] + U[mesh.neighbour()[faceI]]);
        }
    }
    
    Uf = Uf ^ mesh.idir();

    return treconField;
}


tmp<surfaceVectorField> perpVec(const tmp<surfaceScalarField>& tuS)
{
    tmp<surfaceVectorField> tUperp
    (
        TRiSK::perpVec(tuS())
    );
    tuS.clear();
    return tUperp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
