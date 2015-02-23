/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

    For spherical geometry - assumes faces are on a sphere (A spherical face)
    or the area vector points around the sphere (a radial face)

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "VectorSpaceFunctions.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const scalar RS_TOL = 1e-7;

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcSphFaceCentreAndArea
(
    const pointField& p,
    vector& fCtr,
    vector& fArea,
    const bool setCentre
) const
{
    label nPoints = p.size();
    
    // Initial approximate centre and area
    fArea = (p[1]-p[0])^(p[2]-p[0]);
    if (setCentre) // && magSqr(fCtr)<SMALL)
    {
        fCtr = p[0];
        for(int i = 1; i < p.size(); i++) {fCtr += p[i];}
        fCtr /= p.size();
    }

    // determine if this is a spherical or radial face using old vectors
    const scalar cosFaceAngle = (fArea & fCtr) / (mag(fArea)* mag(fCtr));

    const bool faceOnSphere = mag(1 - cosFaceAngle) < 0.5
                           || mag(1 + cosFaceAngle) < 0.5;
        
    // Area and centre for a radial face
    if (!faceOnSphere)
    {
        // Check that face has 4 vertices
        if (nPoints != 4)
        {
            Pout << "faceAngle = " << acos(cosFaceAngle) << " faceAreas = "
                 << fArea << " faceCentres = "
                 << fCtr << endl;
            
            FatalErrorIn("primitiveMesh::calcSphFaceCentreAndArea")
                << " radial faces should have 4 vertices, but face "
                 << "has " << nPoints << exit(FatalError);
        }
        
        // Determine the 4 radii
        scalar a = mag(p[0]);
        scalar b = mag(p[1]);
        scalar c = mag(p[2]);
        scalar d = mag(p[3]);
        
        // First determine the 2 radial directions
        vector rhata = p[0]/a;
        vector rhatb = p[1]/b;
        vector rhatc = p[2]/c;
        vector rhatd = p[3]/d;
            
        // We require b>a and c>d with a and b together and c and d together
        // First check if a and d have the same radial direction
        if (magSqr(rhata - rhatd) <= RS_TOL*magSqr(rhata - rhatc))
        {
            // rotate points around by 1
            scalar dOld = d;
            d = c;
            c = b;
            b = a;
            a = dOld;
            vector rhatdOld = rhatd;
            rhatd = rhatc;
            rhatc = rhatb;
            rhatb = rhata;
            rhata = rhatdOld;
        }
        // Next check if a > b then move points around by 2
        if (a > b)
        {
            scalar dOld = d;
            scalar cOld = c;
            d = b;
            c = a;
            b = dOld;
            a = cOld;
            vector rhatdOld = rhatd;
            vector rhatcOld = rhatc;
            rhatd = rhatb;
            rhatc = rhata;
            rhatb = rhatdOld;
            rhata = rhatcOld;
        }
        const vector rhat0 = rhatc;
        const vector rhat1 = rhata;

        // Check that there are only 2 different rhat
        const scalar distS = magSqr(rhat0 - rhat1);
        if
        (
            magSqr(rhata - rhatb) >= distS*RS_TOL
         || magSqr(rhatc - rhatd) >= distS*RS_TOL
        )
        {
            FatalErrorIn("primitiveMesh::calcSphFaceCentreAndArea")
              << " radial faces should have only 2 radial directions but face "
              << "has\n" << rhata << nl << rhatb << nl << rhatc << nl << rhatd
              << "\ncosFaceAngle = " << cosFaceAngle
              << "\nface points\n" 
              << p[0] << nl << p[1] << nl << p[2] << nl << p[3] << nl
              << "fArea = " << fArea << nl << "fCtr = " << fCtr
              << exit(FatalError);
        }
            
        const vector rhat = setCentre ? unitVector(rhat0 + rhat1)
                                      : unitVector(fCtr);
        fCtr = 0.25*(a + b + c + d)*rhat;
        const vector idir = rhat ^ (rhat1 - rhat0);
        const scalar A = arcLength(rhat0, rhat1)/6.*
        (
            sqr(c) + b*c + sqr(b) - (sqr(a) + a*d + sqr(d))
            //3*(b*c - a*d) - sqr(d-a) + sqr(c-b)
        );
        fArea = A*idir/mag(idir);
        const vector dir = (p[1] - p[0])^(p[2] - p[0]);
        if ((dir & idir) < 0) fArea = -fArea;
    }
    else // face on sphere
    {
        for(int i = 0; i < 2; i++)
        {
            // Use an initial approximate face centre
            point fCentre = fCtr;

            // next find exact centre and the total solid angle, omega
            // don't use solid angle, use projected triangle area
            vector sumAc = vector::zero;
            scalar sumR = 0;
            scalar omega = 0;
            vector sumAn = vector::zero;

            for (label ptI = 0; ptI < nPoints; ptI++)
            {
                const point& nextPoint = p[(ptI + 1) % nPoints];

                vector c = p[ptI] + nextPoint + fCentre;
                scalar r = mag(p[ptI]) + mag(nextPoint);
                vector n = (nextPoint - p[ptI])^(fCentre - p[ptI]);
                scalar a = sphTriSolidAngle(p[ptI], nextPoint, fCentre);
                //scalar a = 0.5*mag(n)/sqr(0.5*r);
                    
                omega += a;
                sumAc += a*c;
                sumR += r*a;
                sumAn += n;
            }
            if (mag(omega) < SMALL)
            {
                FatalErrorIn("calcSphFaceCentreAndArea")
                    << " face " << " has a solid anlge " << omega
                    << exit(FatalError);
            }
                
            sumR /= (2*omega);
            vector rhat = setCentre ? sumAc/mag(sumAc)
                                    : unitVector(fCtr);
            fCtr = rhat*sumR;
            fArea = sqr(sumR)*rhat*omega*sign(sumAn & rhat);
        }
    }
}

void primitiveMesh::calcSphFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas,
    const bool setCentres
) const
{
    // New code for orography
    // Loop through every face
    forAll (faces(), facei)
    {
        const labelList& f = faces()[facei];
        pointField facePoints(f.size());
        forAll(facePoints, i)
        {
            facePoints[i] = p[f[i]];
        }
        
        calcSphFaceCentreAndArea
        (
            facePoints, fCtrs[facei], fAreas[facei], setCentres
        );
    }
}

scalar primitiveMesh::calcSphCellVol
(
    const vectorField& fCtrs,
    const vectorField& fAreas
) const
{
    scalar cellVol = 0;
    
    // approximate cell centre to check face areas for being outward pointing
    point C = vector::zero;
    scalar Cr = 0;
    forAll(fCtrs, facei)
    {
        C += fCtrs[facei];
        Cr += mag(fCtrs[facei]);
    }
    C /= fCtrs.size();
    C = Cr/fCtrs.size()*unitVector(C);

    // Accumulate cell volume
    forAll(fCtrs, facei)
    {
        int outward = sign(fAreas[facei] & (fCtrs[facei] - C));
        scalar SfdotCf = fAreas[facei] & fCtrs[facei];
        cellVol += outward*SfdotCf;
        //Info << "outward = " << outward << " SfdotCf = " << SfdotCf << nl;
    }
    return cellVol/3.; 
}


void primitiveMesh::calcSphCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarList& cellVols,
    const bool setCentres
) const
{
    // New code for with orography
    // Clear the fields for accumulation
    const vectorField oldCellCtrs = cellCtrs;
    cellCtrs = vector::zero;
    cellVols = 0.0;
    
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    
    // Sum of solid angles, alpha, for each cell
    scalarField alphaSum(nCells(), scalar(0));
    
    forAll(own, facei)
    {
        // Calculate volume of sector of the sphere
        scalar SfdotCf = fAreas[facei] & fCtrs[facei];
        // Accumulate onto owner and neighbour cells
        cellVols[own[facei]] += SfdotCf;
        if (facei < nInternalFaces())
        {
            cellVols[nei[facei]] -= SfdotCf;
        }
        
        // Calculate solid angle, alpha
        scalar alpha = mag(SfdotCf)/pow(magSqr(fCtrs[facei]), 1.5);
        
        // Accumulate solid angles and centres for owner and neighbour cells
        cellCtrs[own[facei]] += alpha*fCtrs[facei];
        alphaSum[own[facei]] += alpha;
        if (facei < nInternalFaces())
        {
            cellCtrs[nei[facei]] += alpha*fCtrs[facei];
            alphaSum[nei[facei]] += alpha;
        }
    }
    
    const scalar third = 1./3.;
    forAll(cellVols, i) { cellVols[i] *= third; }
    cellCtrs /= alphaSum;
    
    if (!setCentres)
    {
        cellCtrs = oldCellCtrs/mag(oldCellCtrs)*mag(cellCtrs);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
