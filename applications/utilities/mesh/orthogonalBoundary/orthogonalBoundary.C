/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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

Application
    orthogonalBoundary

Description
    Iteratively moves the points on the boundary in order to make the boundary
    faces orthogonal

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    // New point locations after boundary points have been moved
    IOField<point> newPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh", mesh),
        mesh.points()
    );
    
    // Convergence testing variables
    const label maxIter = 10;
    bool converged = false;
    const scalar moveMax = 1e-6;
    
    for(label iter = 0; iter < maxIter && !converged; iter++)
    {
        scalar maxPointMove = 0;
        // Loop over all boundary faces and move the boundary points
        for(label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
        {
            // Face centre and the projectionPoint (the intersection between
            // face normal and the face, starting from the cell centre)
            const vector& Cf = mesh.faceCentres()[facei];
            const vector& C = mesh.cellCentres()[mesh.faceOwner()[facei]];
            const scalar magSf = mag(mesh.faceAreas()[facei]);
            const vector Sfhat = mesh.faceAreas()[facei]/magSf;
            const vector projPt = C + ((Cf - C) & Sfhat)*Sfhat;
        
            // the face size (used for scaling maxPointMove)
            const scalar faceLength = Foam::sqrt(magSf);
            
            // vector to move the points
            const vector movePoints = projPt - Cf;
            
            // used to test convergenced
            if (mag(movePoints)/faceLength > maxPointMove)
            {
                maxPointMove = mag(movePoints)/faceLength;
            }
            
            // Move the points
            for(label ipp = 0; ipp < mesh.faces()[facei].size(); ipp++)
            {
                label ip = mesh.faces()[facei][ipp];
                newPoints[ip] += movePoints;
            }
        }
        mesh.movePoints(newPoints);
        Info << "maxPointMove = " << maxPointMove << endl;
        converged = (maxPointMove <= moveMax);
    }
    
    newPoints.write();
}
