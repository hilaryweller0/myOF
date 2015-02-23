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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::TRiSKData::calcCircToFluxWeights()
{
    const fvMeshWithDual& mesh = smesh();
    //const fvMeshWithDual& dualMesh = mesh.dualMesh();
    //const polyPatch& bottomPatch = mesh.bottomPatch();
    //const polyPatch& dualPatch = dualMesh.bottomPatch();

    List<FixedList<scalar, 6> >& H = circToFluxOffDiagWeights_;
    List<FixedList<label, 6> >& stencil = circToFluxStencil_;
    labelList& stencilSize = circToFluxStencilSize_;
    surfaceScalarField& Hdiag = circToFluxDiag_;
    
    labelListList& Hp = ddirToFluxOffDiagWeights_;
    labelListList& stencilp = ddirToFluxStencil_;
    
    forAll(Hdiag, faceI)
    {
        if (Htype() == DIAGONAL)
        {
            Hdiag[faceI] = 1;
            H[faceI] = scalar(0);
            stencil[faceI] = label(0);
            stencilSize[faceI] = 0;
        }
        // 3d stencils (using primal only)
        if (Htype() == ASYMMETRIC)
        {
            reconstructH(Hp[faceI], stencilp[faceI], Hdiag[faceI], faceI);
        }
        
        // 2d stencils use the dual grid (vertical faces only)
        const label faced = mesh.dualFaceMap()[faceI];
        if (faced != -1)
        {
            // Determine which H to calculate for this face
        
            // Use Dubos H operator if owner and neibour in dual mesh both quads
            //bool bothQuads=dualMesh.cells()[dualMesh.owner()[faced]].size()==6
            //       && dualMesh.cells()[dualMesh.neighbour()[faced]].size()==6;

            if (Htype() == DUBOS)
            {
                stencilSize[faceI] = 4;
                calcDubosH(H[faceI], stencil[faceI], Hdiag[faceI], faceI);
            }
            else if (Htype() == DUBOSW)
            {
                stencilSize[faceI] = 4;
                calcDubosWH(H[faceI], stencil[faceI], Hdiag[faceI], faceI);
            }
//            else if(Htype() == ASYMMETRICCOMPACT)
//            {
//                stencilSize[faceI] = 4;
//                bool compact = true;
//                reconstructH
//                (
//                    H[faceI], stencil[faceI], Hdiag[faceI], faceI, compact
//                );
//            }
            else if(Htype() == ASYMMETRIC)
            {
                Hdiag[faceI] = mesh.Hdiag()[faceI];
//                stencilSize[faceI] = 6;
//                bool compact = false;
//                reconstructH
//                (
//                    H[faceI], stencil[faceI], Hdiag[faceI], faceI, compact
//                );
            }
        }
    }
    
//    // Make H symmetric if necessary
//    if (Htype() != ASYMMETRIC && Htype() != DUBOS && Htype() != DIAGONAL)
//    {
//        symmetriseH();
//    }
}


void Foam::TRiSKData::calcDubosH
(
    FixedList<scalar,6>& H,
    FixedList<label,6>& stencil,
    scalar& Hdiag,
    const label faceI
)
{
    const fvMeshWithDual& mesh = smesh();
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    //const polyPatch& bottomPatch = mesh.bottomPatch();
    const polyPatch& dualPatch = dualMesh.bottomPatch();
    
    Hdiag = 0;
    
    const label faced = mesh.dualFaceMap()[faceI];

    // The vector for this face
    const vector de = dualMesh.magSf()[faced]*dualMesh.jdir()[faced];

    // consider all faces of the owner and neighbour dual cells
    label celli = dualMesh.owner()[faced];
    bool ownerNeibDone = false;
    label is = 0;
    while(!ownerNeibDone)
    {
        const cell& c = dualMesh.cells()[celli];
        const label fep = (c.size() == 5) ? 6 :(c.size() == 6) ? 4 : 0;
        
        // the vertices at either end of the face edge
        const label ie = dualMesh.faceToPatchEdge()[faced];
        const label iv0 = dualPatch.edges()[ie][0];
        const label iv1 = dualPatch.edges()[ie][1];
        
        // circulate around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label fdi = c[fi];
            // the vertices at either end of the face edge
            const label iee = dualMesh.faceToPatchEdge()[fdi];
            const label iv2 = dualPatch.edges()[iee][0];
            const label iv3 = dualPatch.edges()[iee][1];

            // only use if it is vertical and is adjacent to faced
            if
            (
                dualMesh.dualFaceMap()[fdi] != -1
             && fdi != faced
             && (iv2 == iv0 || iv2 == iv1 || iv3 == iv0 || iv3 == iv1)
            )
            {
                vector dep = dualMesh.magSf()[fdi]*dualMesh.jdir()[fdi];
                scalar denom = fep*mag(de ^ dep);
                
                Hdiag -= (dep & dep)/denom;
                H[is] = (de & dep)*mesh.signMap(faceI)/denom;
                stencil[is] = fdi;
                is++;
            }
        }
        
        if (celli == dualMesh.neighbour()[faced]) ownerNeibDone = true;
        else if (faced < dualMesh.nInternalFaces())
        {
            celli = dualMesh.neighbour()[faced];
        }
        else ownerNeibDone = true;
    }
    if (is != 4)
    {
        FatalErrorIn("dData::calcDubosH")
            << " face " << faceI << " only has " << is
            << " face neighbours in the dual" << exit(FatalError);
    }
    
    // correct H diag for fluxes
//    scalar idotj = mesh.idir()[faceI] & dualMesh.jdir()[faced];
    scalar Sp = mesh.magSf()[faceI];
    scalar Sd = dualMesh.magSf()[faced];
//    Hdiag = -sign(idotj) * Sd/Sp*Hdiag;
    Hdiag *= -Sd/Sp;
}

void Foam::TRiSKData::calcDubosWH
(
    FixedList<scalar,6>& H,
    FixedList<label,6>& stencil,
    scalar& Hdiag,
    const label faceI
)
{
    const fvMeshWithDual& mesh = smesh();
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    //const polyPatch& bottomPatch = mesh.bottomPatch();
    const polyPatch& dualPatch = dualMesh.bottomPatch();
    
    Hdiag = 0;
    
    const label faced = mesh.dualFaceMap()[faceI];

    // The vector for this face
    const vector de = dualMesh.magSf()[faced]*dualMesh.jdir()[faced]
 /(dualMesh.depthf()[faced]*(dualMesh.heightf()[faced]+dualMesh.earthRadius().value()));

    // The edge of the bottom patch associated with this face
    const label ie = dualMesh.faceToPatchEdge()[faced];
    // the vertices at either end of the face edge
    const label iv0 = dualPatch.edges()[ie][0];
    const label iv1 = dualPatch.edges()[ie][1];
    // The area associated with this edge
    const scalar Ae = dualMesh.edgeAreaOwn()[ie] + dualMesh.edgeAreaNei()[ie];

    // consider all faces of the owner and neighbour dual cells
    label celli = dualMesh.owner()[faced];
    label ci = dualMesh.cellToPatchFace()[celli];
    bool ownerNeibDone = false;
    label is = 0;
    // The edge area for the owner/neighbour (owner first)
    scalar Ave = dualMesh.edgeAreaOwn()[ie];
    while(!ownerNeibDone)
    {
        const cell& c = dualMesh.cells()[celli];
        
        // circulate around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label fdi = c[fi];
            // the vertices at either end of the face edge
            const label iee = dualMesh.faceToPatchEdge()[fdi];
            const label iv2 = dualPatch.edges()[iee][0];
            const label iv3 = dualPatch.edges()[iee][1];

            // only use if it is vertical and is adjacent to faced
            if
            (
                dualMesh.dualFaceMap()[fdi] != -1
             && fdi != faced
             && (iv2 == iv0 || iv2 == iv1 || iv3 == iv0 || iv3 == iv1)
            )
            {
                // find the vertex in common
                label ivcommon = iv0;
                if (iv1 == iv2 || iv1 == iv3) ivcommon = iv1;
            
                // Area associated with primal and dual cells
                // 4 points of the kite shape
                point pivc = dualPatch.localPoints()[ivcommon];
                point Cd = dualPatch.faceCentres()[ci];
                point Cf = unitVector(mesh.faceCentres()[faceI]);
                point Cfi = unitVector(dualMesh.faceCentres()[fdi]);
                scalar Aiv = sphTriDistAngle(pivc, Cd, Cf)
                           + sphTriDistAngle(pivc, Cd, Cfi);
                
                vector dep = dualMesh.magSf()[fdi]*dualMesh.jdir()[fdi]
     /(dualMesh.depthf()[fdi]*(dualMesh.heightf()[fdi]+dualMesh.earthRadius().value()));
                scalar denom = magSqr(de ^ dep);
                
                if (faceI == 7386 || faceI == 7389)
                {
                Info << "Ave = " << Ave << " Aiv = " << Aiv << " Ae = " << Ae
                     << " mag(dep) = " << mag(dep) << " magSqr(de ^ dep) = "
                     << magSqr(de ^ dep) << endl;
                }
                Hdiag -= 2*Ave*Aiv/Ae*(dep & dep)/denom;
                H[is] = 2*Ave*Aiv/Ae*(de & dep)*mesh.signMap(faceI)/denom;
                stencil[is] = fdi;
                is++;
            }
        }
        
        if (celli == dualMesh.neighbour()[faced]) ownerNeibDone = true;
        else if (faced < dualMesh.nInternalFaces())
        {
            celli = dualMesh.neighbour()[faced];
            Ave = dualMesh.edgeAreaNei()[ie];
        }
        else ownerNeibDone = true;
    }
    if (is != 4)
    {
        FatalErrorIn("dData::calcDubosH")
            << " face " << faceI << " only has " << is
            << " face neighbours in the dual" << exit(FatalError);
    }
    
    // correct H diag for fluxes
//    scalar idotj = mesh.idir()[faceI] & dualMesh.jdir()[faced];
    scalar Sp = mesh.magSf()[faceI];
    scalar Sd = dualMesh.magSf()[faced];
//    Hdiag = -sign(idotj) * Sd/Sp*Hdiag;
    Hdiag *= -Sd/Sp;
    
    if (faceI == 7386 || faceI == 7389)
    {
        Info << "Face " << faceI << " Hdiag  " << Hdiag << " H = " << H << nl;
    }
}


void Foam::TRiSKData::reconstructH
(
    FixedList<scalar,6>& H,
    FixedList<label,6>& stencil,
    scalar& Hdiag,
    const label faceI,
    const bool compact
)
{
    const fvMeshWithDual& mesh = smesh();
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    //const polyPatch& bottomPatch = mesh.bottomPatch();
    const polyPatch& dualPatch = dualMesh.bottomPatch();

    const label faced = mesh.dualFaceMap()[faceI];

    Hdiag = mag(dualMesh.jdir()[faced] & mesh.idir()[faceI]);
    vector ijdir = mesh.Sf()[faceI] - dualMesh.jdir()[faced]
                  *(dualMesh.jdir()[faced] & mesh.Sf()[faceI]);
    
    // consider all faces of the owner and neighbour dual cells
    // First for calculating tensor A
    symmTensor A(symmTensor::zero);
    
    label celli = dualMesh.owner()[faced];
    bool ownerNeibDone = false;
    label is = 0;
    while(!ownerNeibDone)
    {
        const cell& c = dualMesh.cells()[celli];
        
        // the vertices at either end of the face edge
        const label ie = dualMesh.faceToPatchEdge()[faced];
        const label iv0 = dualPatch.edges()[ie][0];
        const label iv1 = dualPatch.edges()[ie][1];

        // circulate around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label fdi = c[fi];

            // the vertices at either end of the face edge
            const label iee = dualMesh.faceToPatchEdge()[fdi];
            const label iv2 = dualPatch.edges()[iee][0];
            const label iv3 = dualPatch.edges()[iee][1];

            // only use if it is vertical (and is adjacent to faced for compact)
            if
            (
                dualMesh.dualFaceMap()[fdi] != -1
             //&& fdi != faced
             && (!compact ||(iv2 == iv0 || iv2 == iv1 || iv3 == iv0 ||iv3==iv1))
            )
            {
                A += dualMesh.magSf()[fdi]*sqr(dualMesh.jdir()[fdi]);
                is++;
            }
        }
        
        if (celli == dualMesh.neighbour()[faced]) ownerNeibDone = true;
        else if (faced < dualMesh.nInternalFaces())
        {
            celli = dualMesh.neighbour()[faced];
        }
        else ownerNeibDone = true;
    }
    A = inv(A);

    // now calculate the H matrix coefficients for non-quad cells
    Hdiag -= 2*((A & dualMesh.jdir()[faced]) & ijdir)
            *dualMesh.magSf()[faced]/mesh.magSf()[faceI];
    celli = dualMesh.owner()[faced];
    ownerNeibDone = false;
    is = 0;
    while(!ownerNeibDone)
    {
        const cell& c = dualMesh.cells()[celli];
        
        // the vertices at either end of the face edge
        const label ie = dualMesh.faceToPatchEdge()[faced];
        const label iv0 = dualPatch.edges()[ie][0];
        const label iv1 = dualPatch.edges()[ie][1];

        // circulate around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label fdi = c[fi];

            // the vertices at either end of the face edge
            const label iee = dualMesh.faceToPatchEdge()[fdi];
            const label iv2 = dualPatch.edges()[iee][0];
            const label iv3 = dualPatch.edges()[iee][1];

            // only use if it is vertical (and is adjacent to faced)
            if
            (
                (dualMesh.dualFaceMap()[fdi] != -1)
             && (fdi != faced)
             && (!compact ||(iv2 == iv0 || iv2 == iv1 || iv3 == iv0 ||iv3==iv1))
            )
            {
                H[is] = ((A & dualMesh.jdir()[fdi]) & ijdir);
                stencil[is] = fdi;
                is++;
            }
        }
        
        if (celli == dualMesh.neighbour()[faced]) ownerNeibDone = true;
        else if (faced < dualMesh.nInternalFaces())
        {
            celli = dualMesh.neighbour()[faced];
        }
        else ownerNeibDone = true;
    }

    if ((compact && is != 4) || (!compact && (is < 4 || is > 6)))
    {
        FatalErrorIn("TRiSKData::reconstructH")
            << " face " << faceI << " has " << is
            << "elements in the H stencil " << exit(FatalError);
    }
}


void Foam::TRiSKData::reconstructH
(
    labelList& Hp,
    labelList& stencilp,
    scalar& Hdiag,
    const label faceI
)
{
    const fvMeshWithDual& mesh = smesh();
    
    label stencilSize = mesh.cells()[mesh.owner()[faceI]].size()-1;
    if (faceI < mesh.nInternalFaces())
    {
        stencilSize += mesh.cells()[mesh.neighbour()[faceI]].size()-1;
    }
    if (mesh.nGeometricD() == 2) stencilSize -= 4;
    
    Hp.setSize(stencilSize);
    stencilp.setSize(stencilSize);

    Hdiag = mesh.Hdiag()[faceI]; //mesh.ddir()[faceI] & mesh.idir()[faceI];
    vector ijdir = mesh.Sf()[faceI] - mesh.ddir()[faceI]
                  *(mesh.ddir()[faceI] & mesh.Sf()[faceI]);
    
    // consider all faces of the owner and neighbour cells
    // First for calculating tensor A
    symmTensor A = symmTensor::zero;
    
    label celli = mesh.owner()[faceI];
    bool ownerNeibDone = false;
    while(!ownerNeibDone)
    {
        const cell& c = mesh.cells()[celli];
        
        // circulate around the faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label facej = c[fi];
            // Do not include faceI or empty faces
            bool useFace = facej != faceI
                         && facej < mesh.nInternalFaces();
            
            if (!useFace && facej != faceI)
            {
                const label patchID = mesh.boundaryMesh().whichPatch(facej);
                useFace = !isA<emptyPolyPatch>(mesh.boundaryMesh()[patchID]);
            }
            
            if (useFace)
            {
                A += mag(mesh.faceAreas()[facej])
                    *sqr(fieldAccess(mesh.ddir(), facej));
            }
        }
        
        if (celli == mesh.neighbour()[faceI]) ownerNeibDone = true;
        else if (faceI < mesh.nInternalFaces())
        {
            celli = mesh.neighbour()[faceI];
        }
        else ownerNeibDone = true;
    }
    A = inv(A);

    // now calculate the H matrix coefficients
    celli = mesh.owner()[faceI];
    ownerNeibDone = false;
    label is = 0;
    while(!ownerNeibDone)
    {
        const cell& c = mesh.cells()[celli];
        
        // circulate around the faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label facej = c[fi];
            // Do not include faceI or empty faces
            bool useFace = facej != faceI
                         && facej < mesh.nInternalFaces();
            
            if (!useFace && facej != faceI)
            {
                const label patchID = mesh.boundaryMesh().whichPatch(facej);
                useFace = !isA<emptyPolyPatch>(mesh.boundaryMesh()[patchID]);
            }
            
            if (useFace)
            {
                Hp[is] = ((A & mesh.ddir()[facej]) & ijdir);
                stencilp[is] = facej;
                is++;
            }
        }
        
        if (celli == mesh.neighbour()[faceI]) ownerNeibDone = true;
        else if (faceI < mesh.nInternalFaces())
        {
            celli = mesh.neighbour()[faceI];
        }
        else ownerNeibDone = true;
    }

    if (is != stencilSize)
    {
        FatalErrorIn("TRiSKData::reconstructH")
            << " face " << faceI << " has stencil size " << is
            << " should be " << stencilSize << exit(FatalError);
    }
}


void Foam::TRiSKData::symmetriseH()
{
    const fvMeshWithDual& mesh = smesh();
    const fvMeshWithDual& dualMesh = mesh.dualMesh();
    //const polyPatch& bottomPatch = mesh.bottomPatch();
    //const polyPatch& dualPatch = dualMesh.bottomPatch();

    List<FixedList<scalar, 6> >& H = circToFluxOffDiagWeights_;
    List<FixedList<label, 6> >& stencil = circToFluxStencil_;

    Info << "Symmetrising circToFluxOffDiagWeights" << endl;
    List<FixedList<scalar, 6> > Ht(mesh.nFaces());
    forAll(H, faceI)
    {
        // use the dual grid (vertical faces only)
        const label faced = mesh.dualFaceMap()[faceI];
        if (faced != -1)
        {
            // for each H[faceI][is] set j = stencil[faceI][is] and find iss s.t.
            // faceI = stencil[j][iss] so that Ht[faceI][is] = H[j][iss]
            forAll(stencil[faceI], is)
            {
                label jd = stencil[faceI][is];
                label jp = dualMesh.dualFaceMap()[jd];
                label iss = -1;
                for(label i = 0; iss == -1 && i < stencil[jp].size(); i++)
                {
                    if (stencil[jp][i] == faced) iss = i;
                }
                if (iss == -1)
                {
                    FatalErrorIn("dData::calcCircToFluxWeights")
                        << "cannot find " << faced << " in stencil "
                        << stencil[jp] << " for primal face " << faceI
                        << abort(FatalError);
                }
                Ht[faceI][is] = mesh.signMap(jp)*mesh.signMap(faceI)*H[jp][iss];
            }
        }
    }
    // Set H = 0.5*(H + Ht) so that it is symmetric
    forAll(H, faceI)
    {
        if (mesh.dualFaceMap()[faceI] != -1)
        {
            forAll(stencil[faceI],is)
            {
                H[faceI][is] = 0.5*(H[faceI][is] + Ht[faceI][is]);
            }
        }
    }
}
