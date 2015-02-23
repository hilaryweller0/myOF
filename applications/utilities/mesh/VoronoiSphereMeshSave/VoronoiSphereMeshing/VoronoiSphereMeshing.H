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

Class
    VoronoiSphereMeshing

Description
    Conformal-Voronoi automatic mesher for the surface of a sphere
    with grid or read initial points and point position relaxation

SourceFiles
    CGALTriangulation3Ddefs.H
    indexedVertex.H
    indexedCell.H
    VoronoiSphereMeshingI.H
    VoronoiSphereMeshing.C
    VoronoiSphereMeshingIO.C
    VoronoiSphereMeshingAlgorithms.C
    VoronoiSphereMeshingHelpers.C
    VoronoiSphereMeshingNav.C
    
\*---------------------------------------------------------------------------*/

#ifndef VoronoiSphereMeshing_H
#define VoronoiSphereMeshing_H

#define CGAL_INEXACT
#define CGAL_HIERARCHY

#include "CGALTriangulation3Ddefs.H"
#include "meshWithDual.H"
//#include "polyMesh.H"
#include "Switch.H"
#include "InitialPoints.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VoronoiSphereMeshing Declaration
\*---------------------------------------------------------------------------*/

class VoronoiSphereMeshing
:
    public HTriangulation
{
public:

    //- The initial points
    const InitialPoints& initialPoints;
    
    //- The radius of the sphere
    const scalar radius;

    // Control parameters
    #include "VoronoiSphereMeshingControls.H"

private:
    // Private data
       //- list of which locations need smoothing
       boolList meshOk;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        VoronoiSphereMeshing(const VoronoiSphereMeshing&);

        //- Disallow default bitwise assignment
        void operator=(const VoronoiSphereMeshing&);


    // Additional Mesh navigation functions
    #include "VoronoiSphereMeshingNav.H"

public:

    // Constructors

        //- Construct from dictionary
        VoronoiSphereMeshing(const dictionary& dict, const InitialPoints& ips);

    // Destructor

        ~VoronoiSphereMeshing();

    // Member functions
    
        //- Adapt the mesh according to the control dictionary
        void adaptMesh();

        // Conversion functions between point and Point

#           ifdef CGAL_INEXACT
                typedef const Foam::point& pointFromPoint;
                typedef const Point& PointFrompoint;
#           else
                typedef Foam::point pointFromPoint;
                typedef Point PointFrompoint;
#           endif

            inline pointFromPoint topoint(const Point&) const;
            inline PointFrompoint toPoint(const Foam::point&) const;
            inline PointFrompoint toPointWithSize
            (
                const Foam::point&, const scalar
            ) const;

        // Meshing helper functions (in VoronoiSphereMeshingHelpers.C)
        //- The spring force between points x0 and x1
        vector springForce
        (
            const Foam::point& x,
            const Foam::point& y,
            const label np
        ) const;

        //- Find which points need relaxing based on grad(requiredRes)
        void whichPointsOk();

        // Point manipulation

            //- Reindex the points
            void reIndex();

            //- Move the vertex vh to position p
            void movePoint(const Vertex_handle& vh, const Point& p);

        //- Find the intersection point between the facets of the
        //- triangulation and the line between points a and b
        Foam::point intersection
        (
            const Foam::point& a, const Foam::point& b
        ) const;

        // Point insertion and removal

            //- Insert the given points (nearly) on the surface of the sphere.
            void insertPoints(const pointField& points);

            //- Insert the given points close to the given cells (used by
            // addPoints and removeBadShapes and points set to need smoothing
            void insertPoints
            (
                const pointField& points,
                List<Vertex_handle>& vh
            );

            //- remove points in the List (used by removeBadShapes)
             label removePoints(const boolList& toRemove);


        // Algorithms for modifying the mesh
        
            //- Double the resolution globally
            void doubleResolution();

            //- remove points if they are too closely bunched, neighbours
            // set to need smoothing
             label removePoints();

            //- add points if they are too sparse
            // (and set neighbours to need smoothing)
             label addPoints();

            //- Remove squares and shapes with >= 8 sides
            // (neighbours set to need smoothing)
            label removeBadShapes();

            // Make big angles bigger! (thus changing the triangulation
            void enlargeAngles(const label maxEnlargeAnglesIters);

            //- Relax the Tomita (2001) springs
            void TomitaSprings
            (
                const scalar relaxation, const Switch changeTopology = true
            );
            
            void LaplacianSmooth();

            //- Move the points by one Lloyd iteration
            void LloydIteration(const bool realLloyd = true);
            
            //- Move the points to make the mesh more Heikes Randall
            void HRify();

            //- Perform one iteration to move the triangulation towards Pitteway
            // (the return value is the number of non-Pitteway edges)
            label PittewayIteration(const scalar relaxation);

            //- renumber all the points using a CGAL spatial sort
            void spatialSortRenumber();

        // Access

            const pointField dualPoints(const faceList& faces) const;
            
            const pointField triPoints() const;
            
            const faceList dualFaces() const;
            
            const faceList triFaces() const;

//            const scalarList reqRes() const;
            

        // Write

            //- Write the dual or primal mesh
            void write
            (
                const Time& runTime,
                const word& regionName,
                const bool writeDelaunay
            ) const;
            
            //- Write dual patch
            void writePatch(const fileName& fName) const;

            //- Write dual patch
            void writeDualPatch(const fileName& fName) const;

            //- Write the required resolution on the new mesh
//            void writeReqRes(const fileName& fName) const;
};

inline Foam::label sign0(const label i)
    { return i > 0 ? 1 : i == 0 ? 0 : -1; }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VoronoiSphereMeshingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
