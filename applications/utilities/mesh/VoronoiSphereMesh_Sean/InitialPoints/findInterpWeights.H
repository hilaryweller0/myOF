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

Description
    Find the weighting factors for various types of interpolation

Source files
    findInterpWeights.C

\*---------------------------------------------------------------------------*/

#ifndef findInterpWeights_H
#define findInterpWeights_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "point.H"
#include "vectorField.H"
#include "ListOps.H"
#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- The minimum size for an dim dimensional polynomial of order order
inline int minPolySize(const int dim, const int order)
{
    if (order < 0 || order > 2)
    {
        FatalErrorIn("Foam::minSize") << "order must be 0, 1 or 2, not "
            << order << exit(FatalError);
    }

    return order == 0 ? 1 :
           order == 1 ? (dim == 1 ? 2 : dim == 2 ? 3 : dim == 3 ? 4 : 0) :
           order == 2 ? (dim == 1 ? 3 : dim == 2 ? 6 : dim == 3 ? 10 : 0) :
                    0;
}

inline scalar rbf(const scalar r)
{
    return sqrt(sqr(r) + 1);
}

inline void calculateRBFWeights
(
    scalarList& weights,
    const vectorField& points,
    const point& target,
    const label nDims,
    const scalar tol = VSMALL
)
{
    const label N = points.size();
    vectorField pts = points - target;
    scalarList dists2 = magSqr(pts);
    weights.setSize(N);

    // Check for co-incidence with one point
    label closePoint = findMin(dists2);
    // Info << "closePoint = " << closePoint
    //     << " minDist2 = " << dists2[closePoint] << " tol = " << tol << endl;
    if (dists2[closePoint] < tol)
    {
        weights = 0;
        weights[closePoint] = 1;
        // Info << "weights = " << weights << endl;
    }
    else
    {
        // local co-ordinate directions
        label farPoint = findMax(dists2);
        vector idir(pts[farPoint]);
        idir /= mag(idir);
        // next find the furthest point away at 90 degrees to farPoint
        scalar jDist2 = magSqr(pts[0] ^ idir);
        label jPoint = 0;
        for(label j = 1; j < pts.size(); j++)
        {
            if (j != farPoint)
            {
                scalar dist2 = magSqr(pts[j] ^ idir);
                if (dist2 > jDist2)
                {
                    jDist2 = dist2;
                    jPoint = j;
                }
            }
        }
        vector jdir = pts[jPoint] - (pts[jPoint] & idir)*idir;
        jdir /= mag(jdir);

        vector kdir = idir ^ jdir;

        // scale factor for the geometry
        scalar scale = sqrt(max(dists2));

        // Info << "i, j, k = " << idir << " " << jdir << " " <<kdir<<" scale ="
        //     << scale << endl;

        // Re-calculate the list of points in the local co-ordinate directions
        // scaled by the geometry scaling factor.
        forAll(pts, ip)
        {
            const point p = pts[ip];
            pts[ip].x() = (p & idir)/scale;
            pts[ip].y() = nDims >= 2 ? (p & jdir)/scale : 0;
            pts[ip].z() = nDims >= 3 ? (p & kdir)/scale : 0;
        }

        // Info << "pts = " << pts << endl;

        // calculate the matrix of RBF coefficients
        scalarSquareMatrix M(N+1+nDims, N+1+nDims, scalar(0));
        scalarField rbfVals(N + 1 + nDims, scalar(0));
        rbfVals[N] = 1;
        scalar rbf0 = rbf(0);

        for(label i = 0; i < N; i++)
        {
            // Info << "Setting M[" << i << "] of " << N << endl;
            const point& p = pts[i];

            M[i][i] = rbf0;
            M[N][i] = M[i][N] = 1;
            M[N+1][i] = M[i][N+1] = p.x();
            if (nDims >= 2) M[N+2][i] = M[i][N+2] = p.y();
            if (nDims == 3) M[N+3][i] = M[i][N+3] = p.z();

            rbfVals[i] = rbf(mag(p));

            for(label j = i+1; j < N; j++)
            {
                M[i][j] = M[j][i] = rbf(mag(p - pts[j]));
            }
        }

        solve(M, rbfVals);

        // Set the coefficients
        forAll(weights, i)
        {
            weights[i] = rbfVals[i];
        }
    }
}


inline void calculatePolyWeights
(
    scalarList& weights,
    const vectorField& points,
    const point& target,
    const label order,
    const int nDims,
    const label nCentral,
    const scalar centralWeight,
    const scalar tol = VSMALL
)
{
    const int minSize = minPolySize(nDims, order);
    const label N = points.size();
    vectorField pts = points - target;
    scalarList dists2 = magSqr(pts);
    weights.setSize(N);

    // Check for co-incidence with one point
    label closePoint = findMin(dists2);
    if (dists2[closePoint] < tol)
    {
        weights = 0;
        weights[closePoint] = 1;
    }
    else
    {
        // local co-ordinate directions
        label farPoint = findMax(dists2);
        vector idir(pts[farPoint]);
        idir /= mag(idir);
        // next find the furthest point away at 90 degrees to farPoint
        scalar jDist2 = magSqr(pts[0] ^ idir);
        label jPoint = 0;
        for(label j = 1; j < pts.size(); j++)
        {
            if (j != farPoint)
            {
                scalar dist2 = magSqr(pts[j] ^ idir);
                if (dist2 > jDist2)
                {
                    jDist2 = dist2;
                    jPoint = j;
                }
            }
        }
        vector jdir = pts[jPoint] - (pts[jPoint] & idir)*idir;
        jdir /= mag(jdir);

        vector kdir = idir ^ jdir;

        // scale factor for the geometry
        scalar scale = sqrt(max(dists2));

        // Re-calculate the list of points in the local co-ordinate directions
        // scaled by the geometry scaling factor.
        forAll(pts, ip)
        {
            const point p = pts[ip];
            pts[ip].x() = (p & idir)/scale;
            pts[ip].y() = nDims >= 2 ? (p & jdir)/scale : 0;
            pts[ip].z() = nDims >= 3 ? (p & kdir)/scale : 0;
        }

        // weights for different points
        scalarList wts(N, scalar(1));
        for(label iw = 0; iw < nCentral; iw++) wts[iw] = centralWeight;

        // calculate the matrix of the polynomial components
        scalarRectangularMatrix B(N, minSize, scalar(0));

        forAll(pts, ip)
        {
            scalar px = pts[ip].x();
            scalar py = pts[ip].y();
            scalar pz = pts[ip].z();

            label is = 0;
            B[ip][is++] = wts[ip];

            if (order >= 1) B[ip][is++] = wts[ip]*px;
            if (order >= 2) B[ip][is++] = wts[ip]*sqr(px);

            if (nDims >= 2)
            {
                if (order >= 1) B[ip][is++] = wts[ip]*py;
                if (order >= 2)
                {
                    B[ip][is++] = wts[ip]*px*py;
                    B[ip][is++] = wts[ip]*sqr(py);
                }
            }
            if (nDims == 3)
            {
                if (order >= 1) B[ip][is++] = wts[ip]*pz;
                if (order >= 2)
                {
                    B[ip][is++] = wts[ip]*px*pz;
                    B[ip][is++] = wts[ip]*py*pz;
                    B[ip][is++] = wts[ip]*sqr(pz);
                }
            }
        }

        SVD svd(B, SMALL);

        // Set the coefficients
        forAll(weights, i)
        {
            weights[i] = wts[i]*svd.VSinvUt()[0][i];
        }
    }
    if (mag(sum(weights)-1) > 1e-6)
    {
        FatalErrorIn("calculatePolyWeights")
        << " the weights should sum to 1 but weights = " << weights
            << " and 1-sum(weights) = " << 1 - sum(weights)
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
