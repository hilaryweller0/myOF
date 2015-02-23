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
    tanPoints

Description
    Takes the tan of the mesh points in order to create an equal angle cubed
    sphere. Also create cell centres at the bisection of the cell diagonals

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "plane.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    // first move the points
    {
        #include "createMesh.H"

        IOField<point> newPoints
        (
            IOobject("points", mesh.time().constant(), "polyMesh", mesh),
            mesh.points()/mag(mesh.points())
        );
        const scalar piOn4 = 0.25*M_PI;
        
        forAll(newPoints, ip)
        {
            newPoints[ip] = mag(mesh.points()[ip])*point
            (
                Foam::tan(piOn4*newPoints[ip].x()),
                Foam::tan(piOn4*newPoints[ip].y()),
                Foam::tan(piOn4*newPoints[ip].z())
            );
        }
        
        newPoints.write();
    }
    
    // next move the cell and face centres
    {
        #include "createMesh.H"
        // Create cell centres at the bisection of the quad diagonals
        IOField<point> cellCentres
        (
            IOobject("cellCentres", mesh.time().constant(), "polyMesh", mesh),
            mesh.C()
        );
        IOField<point> faceCentres
        (
            IOobject("faceCentres", mesh.time().constant(), "polyMesh", mesh),
            mesh.faceCentres()
        );
        
        // patch of original mesh
        const polyPatch& pPatch(mesh.boundaryMesh()["originalPatch"]);
        
        // Loop over all cells and move the centre to the diagonal bisection of
        // the underlying patch
        forAll(pPatch, pfacei)
        {
            // 4 corners of the patch face
            const point& p0 = pPatch.points()[pPatch[pfacei][0]];
            const point& p1 = pPatch.points()[pPatch[pfacei][1]];
            const point& p2 = pPatch.points()[pPatch[pfacei][2]];
            const point& p3 = pPatch.points()[pPatch[pfacei][3]];
            
            // planes for each diagonal
            plane d0(p0^p2);
            plane d1(p1^p3);
            // intersection of the planes
            plane::ray r = d0.planeIntersect(d1);
            
            // Cell attached to this face
            const label cellI = pPatch.faceCells()[pfacei];
            // Global face number of this face
            const label faceI = pfacei + pPatch.start();
            
            cellCentres[cellI] = unitVector(r.dir())*sign(r.dir()&mesh.C()[cellI])
                                 *mag(cellCentres[cellI]);

            faceCentres[faceI] = unitVector(r.dir())*sign(r.dir()&mesh.C()[cellI])
                                 *mag(faceCentres[pfacei]);

//            Info << "Moving cell centre from " << mesh.C()[cellI] << nl
//                 << "                   to   " << cellCentres[cellI]
//                 << " distance " << mag(mesh.C()[cellI] - cellCentres[cellI])
//                 << endl;
        }
        
        cellCentres.write();
        faceCentres.write();
    }
}
