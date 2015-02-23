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

#include "TRiSKreconstructVec.H"
#include "TRiSKData.H"
#include "meshWithDual.H"
#include "TRiSKperp.H"
#include "fvc.H"
#include "TRiSKconserveInterp.H"
#include "TRiSKFaceToCellMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<surfaceVectorField> divReconstruct(const surfaceScalarField& phi)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(phi.mesh())));
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    
    surfaceScalarField u = phi/mesh.magSf();
    surfaceScalarField v = dualMesh.dualMap
    (
        TRiSK::perp(phi)/dualMesh.magSf()
    );

    surfaceVectorField jdir = dualMesh.dualMap(dualMesh.idir());
    surfaceScalarField idotj = mesh.idir() & jdir;

    tmp<surfaceVectorField> tV
    (
        new surfaceVectorField
        (
            IOobject(phi.name()+"_3", phi.instance(), phi.mesh()),
            mesh.idir()*(u - idotj*v)/(1 - sqr(idotj))
        )
    );
    surfaceVectorField& V = tV();
    
    V += jdir*(v - idotj*u)/(1 - sqr(idotj));
    
    return tV;
}


tmp<surfaceVectorField> divReconstruct(const tmp<surfaceScalarField>& tphi)
{
    tmp<surfaceVectorField> tU
    (
        TRiSK::divReconstruct(tphi())
    );
    tphi.clear();
    return tU;
}


tmp<surfaceVectorField> divReconstruct3d(const surfaceScalarField& phi)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(phi.mesh())));
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    const TRiSKData& TRiSKData = TRiSKData::New(mesh);
    
    // Horizontal contributions to vertical faces (with horizontal normals)
    
    surfaceScalarField u = phi/mesh.magSf();
    surfaceScalarField v = dualMesh.dualMap
    (
        TRiSK::perp(phi)/dualMesh.magSf()
    );
    surfaceVectorField jdir = dualMesh.dualMap(dualMesh.idir());
    surfaceScalarField idotj = mesh.idir() & jdir;

    tmp<surfaceVectorField> tUf
    (
        new surfaceVectorField
        (
            IOobject(phi.name()+"_3", phi.instance(), phi.mesh()),
            mesh.idir()*(u - idotj*v)/(1 - sqr(idotj))
        )
    );
    surfaceVectorField& Uf = tUf();
    
    Uf += jdir*(v - idotj*u)/(1 - sqr(idotj));
    
    // Vertical contributions (w) to vertical faces (with horiztonal normals)
    forAll(Uf, faceI)
    {
        const label stencilSize =TRiSKData.vertToHorizStencilSize()[faceI];
        if (stencilSize != 0)
        {
            const FixedList<label, 4>& stencil
                 = TRiSKData.vertToHorizStencil()[faceI];
            scalar w = 0;
            for(label is = 0; is < stencilSize; is++)
            {
                label faceI = stencil[is];
                w += fieldAccess(u, faceI)
                    *sign(mesh.faceCentres()[faceI] & mesh.faceAreas()[faceI]);
            }
            w /= scalar(stencilSize);
        
            Uf[faceI] += w*mesh.kdir()[faceI];
        }
    }
    
    // Uf for horizontal faces (vertical pointing normals
    volVectorField U("U", faceToCellMap(Uf));
    
    forAll(Uf, faceI)
    {
        if (mesh.faceToPatchEdge()[faceI] == -1)
        {
            Uf[faceI] += U[mesh.owner()[faceI]] + U[mesh.neighbour()[faceI]];
        }
    }
    
    return tUf;
}


tmp<surfaceVectorField> divReconstruct3d(const tmp<surfaceScalarField>& tphi)
{
    tmp<surfaceVectorField> tU
    (
        TRiSK::divReconstruct3d(tphi())
    );
    tphi.clear();
    return tU;
}


tmp<volVectorField> reconstruct
(
    const surfaceScalarField& v,
    const surfaceVectorField& d
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(v.mesh())));

    volSymmTensorField T
    (
        sqr(unitVector(mesh.C()))*fvc::surfaceSum(mesh.magSf())
      + fvc::surfaceSum(mesh.magSf()*sqr(d))
    );

    tmp<volVectorField> tU
    (
        new volVectorField
        (
            IOobject(v.name()+"_3", v.instance(), mesh),
            inv(T) & fvc::surfaceSum(mesh.magSf()*d*v)
        )
    );
    
    return tU;
}


tmp<surfaceVectorField> faceReconstruct
(
    const surfaceScalarField& v,
    const surfaceVectorField& d
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(v.mesh())));

//    midPoint<vector> midInterpv(mesh);
//    midPoint<symmTensor> midInterpt(mesh);

    volSymmTensorField T
    (
        sqr(unitVector(mesh.C()))*fvc::surfaceSum(mesh.magSf())
      + fvc::surfaceSum(mesh.magSf()*sqr(d))
    );

    tmp<surfaceVectorField> tU
    (
        new surfaceVectorField
        (
            IOobject(v.name()+"_3", v.instance(), v.mesh()),
//            inv(midInterpt.interpolate(T))
//          & midInterpv.interpolate(fvc::surfaceSum(mesh.magSf()*d*v))
            inv(conserveInterp(T))
          & conserveInterp(fvc::surfaceSum(mesh.magSf()*d*v))
//            midInterpv.interpolate(inv(T) & fvc::surfaceSum(mesh.magSf()*d*v))
        )
    );
    surfaceVectorField& U = tU();
    
    U += v*d - (U&d)*d;

    return tU;
}


tmp<surfaceScalarField> speedMap
(
    const surfaceScalarField& v,
    const surfaceVectorField& j,
    const surfaceVectorField& i
)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(v.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();
    
    // Check that i is on the primal mesh
    if (i.mesh() != mesh)
    {
        FatalErrorIn("TRiSK::fluxMap") << "i is not on the dual mesh of v"
            << exit(FatalError);
    }
    // Check that v and j are on the same mesh
    if (v.mesh() != j.mesh())
    {
        FatalErrorIn("TRiSK::fluxMap") << "v and j are not on the same mesh"
            << exit(FatalError);
    }

    tmp<surfaceScalarField> tu
    (
        new surfaceScalarField
        (
            IOobject(v.name()+"_i", v.instance(), mesh),
            (dualMesh.dualMap(faceReconstruct(v,j)) & i)
        )
    );

    return tu;
}

tmp<surfaceScalarField> fluxMap
(
    const surfaceScalarField& vS,
    const surfaceVectorField& j,
    const surfaceVectorField& i
)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();

    return mesh.magSf()*speedMap(vS/dualMesh.magSf(), j, i);
}


tmp<surfaceScalarField> circToFlux(const surfaceScalarField& vS)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();
    const TRiSKData& TRiSKData = TRiSKData::New(mesh);

    if (TRiSKData.Htype() == TRiSKData::DIAGONAL)
    {
        return dualMesh.dualFluxMap(vS);
    }

    return TRiSKData.circToFluxDiag()*dualMesh.dualFluxMap(vS)
          + circToFluxOffDiag(vS);
}

tmp<surfaceScalarField> circToFluxOffDiag(const surfaceScalarField& vS)
{
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const fvMeshWithDual& mesh = dualMesh.dualMesh();
    const TRiSKData& TRiSKData = TRiSKData::New(mesh);
    
    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject(vS.name(), vS.instance(), mesh),
            mesh,
            dimensionedScalar("", vS.dimensions(), scalar(0))
        )
    );
    surfaceScalarField& uS = tuS();
    
    if (TRiSKData.Htype() == TRiSKData::DUBOS)
    {
        forAll(uS, faceI)
        {
            const label faced = mesh.dualFaceMap()[faceI];
            if (faced != -1)
            {
                const FixedList<scalar, 6>& weights
                     = TRiSKData.circToFluxOffDiagWeights()[faceI];
                const FixedList<label, 6>& stencil
                     = TRiSKData.circToFluxStencil()[faceI];
                for(label is = 0; is < TRiSKData.circToFluxStencilSize()[faceI]; is++)
                {
                    uS[faceI] += weights[is]*vS[stencil[is]];
                }
            }
        }
    }
    else if (TRiSKData.Htype() == TRiSKData::ASYMMETRIC)
    {
//        surfaceVectorField Uf
//        (
//            inv(fvc::interpolate(fvc::surfaceSum(dualMesh.magSf()*sqr(dualMesh.jdir())), "H"))
//          & fvc::interpolate(fvc::surfaceSum(dualMesh.jdir()*vS), "H")
//        );

        volSymmTensorField T
        (
            fvc::surfaceSum(dualMesh.magSf()*sqr(dualMesh.jdir()))
//          fvc::surfaceSum(sqr(dualMesh.magSf())*sqr(dualMesh.jdir()))
        );

        volSymmTensorField rs = sqr(dualMesh.C()); // /mag(dualMesh.C());
        T.internalField() += rs.internalField();

        surfaceVectorField Uf
        (
            fvc::interpolate
            (
                inv(T)
//            & fvc::surfaceSum(dualMesh.jdir()*vS/dualMesh.deltaCoeffs()),      
//            & fvc::surfaceSum(dualMesh.jdir()*vS*dualMesh.magSf()),      
              & fvc::surfaceSum(dualMesh.jdir()*vS),      
                "H"
            )
        );

        // Take dot product in Sf dir and correct and remove ddir component
        uS = dualMesh.dualMap(Uf)
           & (mesh.Sf() - mesh.ddir()*mesh.Hdiag()*mesh.magSf());

//        // Take dot product in Sf dir and correct and remove diag component
//        uS = (dualMesh.dualMap(Uf) & mesh.Sf())
//           - TRiSKData.circToFluxDiag()*dualMesh.dualFluxMap(vS);
    }
    
    return tuS;
}

tmp<surfaceScalarField> ddirToFlux(const surfaceScalarField& vS)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const TRiSKData& TRiSKData = TRiSKData::New(mesh);

    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject(vS.name(), vS.instance(), mesh),
            vS
        )
    );
    surfaceScalarField& uS = tuS();

    if (TRiSKData.Htype() != TRiSKData::DIAGONAL)
    {
        uS = mesh.Hdiag()*vS + ddirToFluxOffDiag(vS);
    }

    return tuS;
}


tmp<surfaceScalarField> ddirToFluxOffDiag(const surfaceScalarField& vS)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    const TRiSKData& TRiSKData = TRiSKData::New(mesh);
    
    tmp<surfaceScalarField> tuS
    (
        new surfaceScalarField
        (
            IOobject(vS.name(), vS.instance(), mesh),
            mesh,
            dimensionedScalar("", vS.dimensions(), scalar(0))
        )
    );
    surfaceScalarField& uS = tuS();
    
    if (TRiSKData.Htype() != TRiSKData::DIAGONAL)
    {
        forAll(uS, faceI)
        {
            const labelList& weights
                = TRiSKData.ddirToFluxOffDiagWeights()[faceI];
            const labelList& stencil
                = TRiSKData.ddirToFluxStencil()[faceI];
            for(label is = 0; is < stencil.size(); is++)
            {
                uS[faceI] += weights[is]*vS[stencil[is]];
            }
        }
    }
    
    return tuS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
