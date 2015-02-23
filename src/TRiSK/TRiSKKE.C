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

#include "TRiSKKE.H"
#include "TRiSKData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> KE
(
    const surfaceScalarField& uS,
    const surfaceScalarField& vS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
    if (dualMesh == mesh.dualMesh())
    {
        return KEprimalDual(uS, vS);
    }
    else
    {
        return KEprimal(uS, vS);
    }
}


tmp<volScalarField> KEprimalDual
(
    const surfaceScalarField& uS,
    const surfaceScalarField& vS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));
    const fvMeshWithDual& dualMesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));

    tmp<volScalarField> tKE
    (
        new volScalarField
        (
            IOobject("KE("+uS.name()+")", uS.instance(), mesh),
            mesh,
            dimensionedScalar("KE", dimensionSet(0,2,-2,0,0), scalar(0))
        )
    );
    volScalarField& KE = tKE();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(uS, faceI)
    {
        label faced = mesh.dualFaceMap()[faceI];
        
//        scalar k = -mesh.signMap(faceI)*uS[faceI]/mesh.magSf()[faceI]
//                  * vS[faced]/dualMesh.magSf()[faced];
//        KE[owner[faceI]] += k*mesh.faceVolOwn()[faceI];
//        KE[neighbour[faceI]] += k*mesh.faceVolNei()[faceI];

        scalar k = -0.5*mesh.signMap(faceI)*uS[faceI]*vS[faced]
                    /dualMesh.depthf()[faced];
//        scalar lambda = mesh.surfaceInterpolation::weights()[faceI];
//        KE[owner[faceI]]    += k*(1-lambda);
//        KE[neighbour[faceI]]+= k*lambda;
      KE[owner[faceI]]    += k*mesh.faceVolOwn()[faceI]/mesh.faceVol()[faceI];
      KE[neighbour[faceI]]+= k*mesh.faceVolNei()[faceI]/mesh.faceVol()[faceI];
    }
    
//    forAll(mesh.boundary(), patchi)
//    {
//        const labelUList& pFaceCells =
//            mesh.boundary()[patchi].faceCells();

//        const fvsPatchField<scalar>& puS = uS.boundaryField()[patchi];
//        const fvsPatchField<scalar>& pv = v.boundaryField()[patchi];

//        forAll(mesh.boundary()[patchi], facei)
//        {
//            KE[pFaceCells[facei]] += -= 0.5*puS[facei]*pv[?
//        }
//    }

    KE.internalField() /= mesh.V();

    return tKE;
}


tmp<volScalarField> KEprimal
(
    const surfaceScalarField& uS,
    const surfaceScalarField& vS
)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));
    if (mesh != vS.mesh())
    {
        FatalErrorIn("TRiSK::KEprimal")
             << " called for fields on different meshes"
             << exit(FatalError);
    }

    tmp<volScalarField> tKE
    (
        new volScalarField
        (
            IOobject("KE("+uS.name()+")", uS.instance(), mesh),
            mesh,
            dimensionedScalar("KE", dimensionSet(0,2,-2,0,0), scalar(0))
        )
    );
    volScalarField& KE = tKE();
    
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(uS, faceI)
    {
        scalar k = uS[faceI]*vS[faceI]/sqr(mesh.magSf()[faceI]);
        KE[owner[faceI]] += k*mesh.faceVolOwn()[faceI];
        KE[neighbour[faceI]] += k*mesh.faceVolNei()[faceI];
    }

    KE.internalField() /= mesh.V();

    return tKE;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
