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

#include "TRiSKcurl.H"
//#include "TRiSKperp.H"
#include "TRiSKData.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TRiSK
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> curl(const surfaceScalarField& vS)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(vS.mesh())));
//    const fvMeshWithDual& dualMesh = mesh.dualMesh();

    tmp<volScalarField> tcurl
    (
        new volScalarField
        (
            IOobject("curl("+vS.name()+")", vS.instance(), mesh),
            mesh,
            dimensionedScalar("curl", vS.dimensions()/dimVol, scalar(0))
        )
    );
    volScalarField& c = tcurl();
    
    // loop over all the faces and accumulate the contribution to the curl for
    // the vertical faces
    forAll(vS, faceI)
    {
//        const label faced = mesh.dualFaceMap()[faceI];
//        if (faced != -1 && faced < dualMesh.nInternalFaces())
//        {
//            const label ownd = dualMesh.owner()[faced];
//            const label neid = dualMesh.neighbour()[faced];
            const label own = mesh.owner()[faceI];
            const label nei = mesh.neighbour()[faceI];
            
            // contribution to the curl for each dual cell
            const scalar ci = vS[faceI];//*mesh.magSf()[faceI];
    
            c[own] += ci;
            c[nei] -= ci;
//        }
    }
    c.internalField() /= mesh.V();
    
    return tcurl;
}


tmp<volScalarField> curl(const tmp<surfaceScalarField>& tvS)
{
    tmp<volScalarField> tcurl
    (
        TRiSK::curl(tvS())
    );
    tvS.clear();
    return tcurl;
}


tmp<surfaceVectorField> curl3d(const surfaceScalarField& uS)
{
    const fvMeshWithDual& mesh
         = *(dynamic_cast<const fvMeshWithDual*>(&(uS.mesh())));
    //const fvMeshWithDual& dualMesh = mesh.dualMesh();

    const surfaceScalarField u = uS/mesh.magSf();

    //// Initialise curlf to be the vertical component using the dual
    // Initialise curlf to be zero
    tmp<surfaceVectorField> tcurlf
    (
        new surfaceVectorField
        (
            IOobject("curl("+uS.name()+")", uS.instance(), mesh),
            mesh,
            dimensionedVector("curl", uS.dimensions()/dimVol, vector::zero)
        )
    );
    surfaceVectorField& curlf = tcurlf();
    
    // Curl along all edges 
    scalarList curlEdge(mesh.nEdges(), scalar(0));
    forAll(curlEdge, ie)
    {
        // Edge centre
        const point Ce = unitVector
        (
            mesh.points()[mesh.edges()[ie][0]]
          + mesh.points()[mesh.edges()[ie][1]]
        )*0.5*
        (
            mag(mesh.points()[mesh.edges()[ie][0]])
          + mag(mesh.points()[mesh.edges()[ie][1]])
        );
        // vector of the edge
        const vector e = mesh.points()[mesh.edges()[ie][1]]
                       - mesh.points()[mesh.edges()[ie][0]];

        //bool boundaryEdge = false;

        const labelList& edgeFaces = mesh.edgeFaces()[ie];
        scalar Ae = 0;
//        Info << "Edge " << ie << nl;
        forAll(edgeFaces, fi)
        {
            const label faceI = edgeFaces[fi];
            const point& Co = mesh.C()[mesh.owner()[faceI]];
            const point& Cn = faceI < mesh.nInternalFaces() ?
                              mesh.C()[mesh.neighbour()[faceI]] :
                              mesh.intersections()[faceI];
            const vector delta = Co - Cn;
            // radial vector from edge to delta vector centre
            const vector r = mesh.intersections()[faceI] - Ce;
            //need to get the sign right in comparison to the edge direction
            vector A = r ^ delta;
            //curlEdge[ie] += sign(A & e)*u[faceI]*sphDist(Cn,Co);
            curlEdge[ie] -= sign(A & e)*u[faceI]*mag(delta);
            Ae += mag(A);
            if (faceI >= mesh.nInternalFaces())
            {
                fi = edgeFaces.size();
                curlEdge[ie] = 0;
                //boundaryEdge = true;
            }
        }
        curlEdge[ie] /= 0.5*Ae;
    }
    
    // Reconstruct face values of curl from edge values
    forAll(curlf, faceI)
    {
        symmTensor W = sqr(mesh.Sf()[faceI])/mesh.magSf()[faceI];
        const labelList& faceEdges = mesh.faceEdges()[faceI];
        forAll(faceEdges, ie)
        {
            label edgeI = faceEdges[ie];
            const vector e = mesh.points()[mesh.edges()[edgeI][1]]
                           - mesh.points()[mesh.edges()[edgeI][0]];
            W += sqr(e)/mag(e);
            curlf[faceI] += e*curlEdge[edgeI];
        }
        curlf[faceI] = inv(W) & curlf[faceI];
    }
    
//    // Replace vertical component on vertical faces with 2d version
//    volScalarField curlz("curl", curl(mesh.dualFluxMap(uS)));
//    surfaceScalarField curlzf = dualMesh.dualMap(fvc::interpolate(curlz));
//    forAll(curlf, faceI)
//    {
//        const label faced = mesh.dualFaceMap()[faceI];
//        if (faced != -1 && faced < dualMesh.nInternalFaces())
//        {
//            curlf[faceI] += mesh.rHatf()[faceI]*
//            (
//                curlzf[faceI] - (curlf[faceI] & mesh.rHatf()[faceI])
//            );
//        }
//    }
    
    return tcurlf;
}


tmp<surfaceVectorField> curl3d(const tmp<surfaceScalarField>& tuS)
{
    tmp<surfaceVectorField> tcurl
    (
        TRiSK::curl3d(tuS())
    );
    tuS.clear();
    return tcurl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TRiSK

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
