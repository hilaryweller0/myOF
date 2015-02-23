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

\*----------------------------------------------------------------------------*/

#include "CV2D.H"
#include "Random.H"
#include "transform.H"
#include "IFstream.H"
#include "uint.H"
#include "ulong.H"
#include "Pair.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV2D::insertBoundingBox()
{
    Info<< "insertBoundingBox: creating bounding mesh" << endl;
    scalar bigSpan = 10*tols_.span;
    insertPoint(point2D(-bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(-bigSpan, bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(bigSpan, bigSpan), Vb::FAR_POINT);
}


void Foam::CV2D::fast_restore_Delaunay(Vertex_handle vh)
{
    int i;
    Face_handle f = vh->face(), next, start(f);

    do
    {
        i=f->index(vh);
        if (!is_infinite(f))
        {
            if (!internal_flip(f, cw(i))) external_flip(f, i);
            if (f->neighbor(i) == start) start = f;
        }
        f = f->neighbor(cw(i));
    } while (f != start);
}

void Foam::CV2D::external_flip(Face_handle& f, int i)
{
    Face_handle n = f->neighbor(i);

    if
    (
        CGAL::ON_POSITIVE_SIDE 
     != side_of_oriented_circle(n, f->vertex(i)->point())
    ) return;

    flip(f, i);
    i = n->index(f->vertex(i));
    external_flip(n, i);
}

bool Foam::CV2D::internal_flip(Face_handle& f, int i)
{
    Face_handle n = f->neighbor(i);

    if
    (
        CGAL::ON_POSITIVE_SIDE
     != side_of_oriented_circle(n, f->vertex(i)->point())
    )
    {
        return false;
    }

    flip(f, i);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV2D::CV2D
(
    const dictionary& controlDict,
    const querySurface& qSurf
)
:
    HTriangulation(),
    qSurf_(qSurf),
    controls_(controlDict),
    tols_(controlDict, controls_.minCellSize, qSurf.bb()),
    z_((1.0/3.0)*(qSurf_.bb().min().z() + qSurf_.bb().max().z())),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    startOfBoundaryConformPointPairs_(0)
{
    insertBoundingBox();
    insertFeaturePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CV2D::~CV2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV2D::insertPoints
(
    const point2DField& points,
    const scalar nearness
)
{
    Info<< "insertInitialPoints(const point2DField& points): ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    // Add the points and index them
    forAll(points, i)
    {
        const point2D& p = points[i];

        if (qSurf_.wellInside(toPoint3D(p), nearness))
        {
            insert(toPoint(p))->index() = nVert++;
        }
        else
        {
            Warning 
                << "Rejecting point " << p << " outside surface" << endl;
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (controls_.writeInitialTriangulation)
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
        writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV2D::insertPoints(const fileName& pointFileName)
{
    IFstream pointsFile(pointFileName);

    if (pointsFile.good())
    {
        insertPoints(point2DField(pointsFile), 0.5*controls_.minCellSize2);
    }
    else
    {
        FatalErrorIn("insertInitialPoints")
            << "Could not open pointsFile " << pointFileName
            << exit(FatalError);
    }
}


void Foam::CV2D::insertGrid()
{
    Info<< "insertInitialGrid: ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    scalar x0 = qSurf_.bb().min().x();
    scalar xR = qSurf_.bb().max().x() - x0;
    int ni = int(xR/controls_.minCellSize) + 1;
    scalar deltax = xR/ni;

    scalar y0 = qSurf_.bb().min().y();
    scalar yR = qSurf_.bb().max().y() - y0;
    int nj = int(yR/controls_.minCellSize) + 1;
    scalar deltay = yR/nj;

    Random rndGen(1321);
    scalar pert = controls_.randomPurturbation*min(deltax, deltay);

    for (int i=0; i<ni; i++)
    {
        for (int j=0; j<nj; j++)
        {
            point p(x0 + i*deltax, y0 + j*deltay, 0);

            if (controls_.randomiseInitialGrid)
            {
                p.x() += pert*(rndGen.scalar01() - 0.5);
                p.y() += pert*(rndGen.scalar01() - 0.5);
            }

            if (qSurf_.wellInside(p, 0.5*controls_.minCellSize2))
            {
                insert(Point(p.x(), p.y()))->index() = nVert++;
            }
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (controls_.writeInitialTriangulation)
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
        writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV2D::insertSurfacePointPairs()
{
    startOfSurfacePointPairs_ = number_of_vertices();

    if (controls_.insertSurfaceNearestPointPairs)
    {
        insertSurfaceNearestPointPairs();
    }

    if (controls_.writeNearestTriangulation)
    {
        writeFaces("near_allFaces.obj", false);
        writeFaces("near_faces.obj", true);
        writeTriangles("near_triangles.obj", true);
    }

    // Insertion of point-pais for near-points may cause protrusions
    // so insertBoundaryConformPointPairs must be executed last
    if (controls_.insertSurfaceNearPointPairs)
    {
        insertSurfaceNearPointPairs();
    }

    startOfBoundaryConformPointPairs_ = number_of_vertices();
}


void Foam::CV2D::boundaryConform()
{
    if (!controls_.insertSurfaceNearestPointPairs)
    {
        markNearBoundaryPoints();
    }

    // Mark all the faces as SAVE_CHANGED
    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        fit++
    )
    {
        fit->faceIndex() = Fb::SAVE_CHANGED;
    }

    for (label iter=1; iter<=controls_.maxBoundaryConformingIter; iter++)
    {
        label nIntersections = insertBoundaryConformPointPairs
        (
            "surfaceIntersections_" + Foam::name(iter) + ".obj"
        );

        if (nIntersections == 0)
        {
            break;
        }
        else
        {
            Info<< "BC iteration " << iter << ": "
                << nIntersections << " point-pairs inserted" << endl;
        }

        // Any faces changed by insertBoundaryConformPointPairs will now
        // be marked CHANGED, mark those as SAVE_CHANGED and those that
        // remained SAVE_CHANGED as UNCHANGED
        for
        (
            Triangulation::Finite_faces_iterator fit = finite_faces_begin();
            fit != finite_faces_end();
            fit++
        )
        {
            if (fit->faceIndex() == Fb::SAVE_CHANGED)
            {
                fit->faceIndex() = Fb::UNCHANGED;
            }
            else if (fit->faceIndex() == Fb::CHANGED)
            {
                fit->faceIndex() = Fb::SAVE_CHANGED;
            }
        }
    }

    Info<< nl;
}


void Foam::CV2D::removeSurfacePointPairs()
{
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            Finite_vertices_iterator vir = vit;
            vit--;
            remove(vir);
        }
    }
}


void Foam::CV2D::newPoints(const scalar relaxation)
{
    Info<< "newPointsFromVertices: ";

    const vectorField& faceNormals = qSurf_.faceNormals();

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() && !vit->nearBoundary())
        {
            // Current dual-cell defining vertex ("centre")
            point2DFromPoint defVert0 = toPoint2D(vit->point());

            Triangulation::Edge_circulator ec = incident_edges(vit);
            Triangulation::Edge_circulator ecStart = ec;

            // Circulate around the edges to find the first which is not
            // infinite
            do
            {
                if (!is_infinite(ec)) break;
            } while (++ec != ecStart);

            // Store the start-end of the first non-infinte edge
            point2D de0 = toPoint2D(circumcenter(ec->first));

            // Keep track of the maximum edge length^2
            scalar maxEdgeLen2 = 0.0;

            // Keep track of the index of the longest edge
            label edgecd0i = -1;

            // Edge counter
            label edgei = 0;

            do
            {
                if (!is_infinite(ec))
                {
                    // Get the end of the current edge
                    point2D de1 = toPoint2D
                    (
                        circumcenter(ec->first->neighbor(ec->second))
                    );

                    // Store the current edge vector
                    edges[edgei] = de1 - de0;

                    // Store the edge mid-point in the vertices array
                    vertices[edgei] = 0.5*(de1 + de0);

                    // Move the current edge end into the edge start for the
                    // next iteration
                    de0 = de1;

                    // Keep track of the longest edge

                    scalar edgeLen2 = magSqr(edges[edgei]);

                    if (edgeLen2 > maxEdgeLen2)
                    {
                        maxEdgeLen2 = edgeLen2;
                        edgecd0i = edgei;
                    }

                    edgei++;
                }
            } while (++ec != ecStart);

            // Initialise cd0 such that the mesh will align
            // in in the x-y directions
            vector2D cd0(1, 0);

            if (controls_.relaxOrientation)
            {
                // Get the longest edge from the array and use as the primary
                // direction of the coordinate system of the "square" cell
                cd0 = edges[edgecd0i];
            }

            if (controls_.nearWallAlignedDist > 0)
            {
                pointIndexHit pHit = qSurf_.tree().findNearest
                (
                    toPoint3D(defVert0),
                    controls_.nearWallAlignedDist2
                );

                if (pHit.hit())
                {
                    cd0 = toPoint2D(faceNormals[pHit.index()]);
                }
            }

            // Rotate by 45deg needed to create an averaging procedure which
            // encourages the cells to be square
            cd0 = vector2D(cd0.x() + cd0.y(), cd0.y() - cd0.x());

            // Normalise the primary coordinate direction
            cd0 /= mag(cd0)+SMALL;

            // Calculate the orthogonal coordinate direction
            vector2D cd1(-cd0.y(), cd0.x());


            // Restart the circulator
            ec = ecStart;

            // ... and the counter
            edgei = 0;

            // Initialiose the displacement for the centre and sum-weights
            vector2D disp = vector2D::zero;
            scalar sumw = 0;

            do
            {
                if (!is_infinite(ec))
                {
                    // Pick up the current edge
                    const vector2D& ei = edges[edgei];

                    // Calculate the centre to edge-centre vector
                    vector2D deltai = vertices[edgei] - defVert0;

                    // Set the weight for this edge contribution
                    scalar w = controls_.squares ?
                        magSqr(deltai.x()*ei.y() - deltai.y()*ei.x())
                      : 1; //magSqr(ei)*mag(deltai);

                    /*
                    scalar aspect = 1;

                    if (nearWall)
                    {
                        vector2D ei2 = ei;
                        ei2.x() *= 2;
                        // sqr -> aspect ratio of 2
                        aspect = magSqr(ei2);
                    }
                    */

                    //scalar w = 
                    //    //aspect*
                    //    (1 + 0.5*(grade&deltai)/mag(deltai))
                    //   *magSqr(deltai.x()*ei.y() - deltai.y()*ei.x());

                    if (controls_.squares)
                    {
                        // Use the following for an ~square mesh
                        // Find the coordinate contributions for this edge delta
                        scalar cd0deltai = cd0 & deltai;
                        scalar cd1deltai = cd1 & deltai;

                        // Create a "square" displacement
                        if (mag(cd0deltai) > mag(cd1deltai))
                        {
                            disp += (w*cd0deltai)*cd0;
                        }
                        else
                        {
                            disp += (w*cd1deltai)*cd1;
                        }
                    }
                    else
                    {
                        // Use this for a hexagon/pentagon mesh
                        disp += w*deltai;
                    }

                    // Sum the weights
                    sumw += w;
                }
                else
                {
                    FatalErrorIn("CV2D::newPoints() const")
                        << "Infinite triangle found in internal mesh"
                        << exit(FatalError);
                }

                edgei++;

            } while (++ec != ecStart);

            // Calculate the average displacement
            disp /= sumw;

            // Move the point by a fraction of the average displacement
            movePoint(vit, defVert0 + relaxation*disp);
        }
    }

    Info<< endl;
}


void Foam::CV2D::moveInternalPoints(const point2DField& newPoints)
{
    label pointI = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            movePoint(vit, newPoints[pointI++]);
        }
    }
}


Foam::label Foam::CV2D::PittewayIteration(const scalar relax)
{
    Info << "PittewayIteration(relaxation = " << relax << ")" << endl;

    label nNonPittEdges = 0;

    // Lists of new positions for the vertices and whether each vertex should move
    vector2DField newPositions(label(number_of_vertices()-1), vector2D::zero);
    boolList moveVertex(label(number_of_vertices()), false);

    // Loop through edges
    for
    (
        Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        // Find the two vertices connected by this edge and their indices
        const Vertex_handle vi = eit->first->vertex(cw(eit->second));
        const Vertex_handle vj = eit->first->vertex(ccw(eit->second));
        const label ii = vi->index();
        const label ij = vj->index();
        
        // Only operate on these vertices if they haven't already moved
        // and if the are internal
        if
        (
            !moveVertex[ii] && !moveVertex[ij]
         && vi->internalOrBoundaryPoint() && vj->internalOrBoundaryPoint()
        )
        {
            // The points at these vertices
            const point2D pi = toPoint2D(vi->point());
            const point2D pj = toPoint2D(vj->point());
            // midpoint of the edge
            const point2D pijmid = 0.5*(pi + pj);

            // Find the two outlying points from the edge
            // First define the edge in the opposite direction
            Face_handle fNeib;
            int iNeib;
            bool backEdgeExists = is_edge(vj, vi, fNeib, iNeib);
            if (!backEdgeExists)
            {
                FatalErrorIn("PittewayIteration")
                    << " no backward edge between vertices at " << pi
                    << " and " << pj << exit(FatalError);
            }
            Edge backEdge(fNeib,iNeib);
            Foam::Pair<point2D> pOut
            (
                toPoint2D(eit->first->vertex(eit->second)->point()),
                toPoint2D(fNeib->vertex(iNeib)->point())
            );
        
            // The circumcentres of the facets of the cells
            Foam::Pair<point2D> circums;
            circums[0] = toPoint2D(circumcenter(eit->first));
            //circums[1] = circumcenter(eit->first->neighbor(eit->second));
            circums[1] = toPoint2D(circumcenter(fNeib));
            
            // Only operate on this edge if one of the circumcentres is the
            // wrong side of this edge
            bool OneOnZeroSide =((circums[1]-pijmid) & (pOut[0]-pijmid))
                               > SMALL;
            bool ZeroOnOneSide =((circums[0]-pijmid) & (pOut[1]-pijmid))
                               > SMALL;
            
            if (ZeroOnOneSide)
            {
                pOut = pOut.reversePair();
                circums = circums.reversePair();
                OneOnZeroSide = true;
            }
            if (OneOnZeroSide)
            {
                nNonPittEdges += 1;
                moveVertex[ii] = true;
                moveVertex[ij] = true;
                
                // Define the vector between the circumcentres
                const vector2D Ve = circums[0] - circums[1];
                // and a point along the line closer to circums[1]
                const point2D Vme = (1-SMALL)*circums[1] + SMALL*circums[0];
            
                // Find the fraction to move pi and pj towards pOut
                const scalar alphai = max
                (
                    min(relax*((Vme-pi)&Ve)/((pOut[0]-pi)&Ve), 0.5),
                    0
                );
                const scalar alphaj = max
                (
                    min(relax*((Vme-pj)&Ve)/((pOut[0]-pj)&Ve), 0.5),
                    0
                );
                
                // Calculate the new positions
                newPositions[ii] = (1 - alphai)*pi + alphai*pOut[0];
                newPositions[ij] = (1 - alphaj)*pj + alphaj*pOut[0];
                
                // if vi or vj are boundary points then move the other more
                if (vi->ppMaster())
                {
                    // find the intersections between lines pj->pOut[0] and
                    // pi->circums[1]
                    const line<point2D,const point2D&> pjpO(pj,pOut[0]);
                    const line<point2D,const point2D&> pic1(pi, Vme);
                    point2D intersect1;
                    pjpO.nearestDist(pic1, newPositions[ij], intersect1);
                    newPositions[ij] = relax*newPositions[ij] + (1-relax)*pj;
                }
                else if (vj->ppMaster())
                {
                    // find the intersections between lines pi->pOut[0] and
                    // pj->circums[1]
                    const line<point2D,const point2D&> pipO(pi,pOut[0]);
                    const line<point2D,const point2D&> pjc1(pj, Vme);
                    point2D intersect0;
                    point2D intersect1;
                    pipO.nearestDist(pjc1, newPositions[ii], intersect1);
                    newPositions[ii] = relax*newPositions[ii] + (1-relax)*pi;
                }
        
//                Info << "Moving " << pi << " to " << newPositions[ii] << nl
//                     << "  and  " << pj << " to " << newPositions[ij] << nl
//                     << "circums = " << circums << nl
//                     << " pOut  = " << pOut << endl;
            }
        }
    }
    
    // Move only the internal vertices
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        const label iv = vih->index();
        if (moveVertex[iv] && vih->internalPoint())
        {
            movePoint(vih, newPositions[iv]);
        }
    }
    
    Info << "    modified " << nNonPittEdges << " non-Pitteway edges" << endl;

    return nNonPittEdges;

}


void Foam::CV2D::write() const
{
    if (controls_.writeFinalTriangulation)
    {
        writePoints("internalPoints.obj", true);
        writeFaces("allFaces.obj", false);
        writeFaces("faces.obj", true);
        writeTriangles("allTriangles.obj", false);
        writeTriangles("triangles.obj", true);
        writePatch("patch.pch");
    }
}


// ************************************************************************* //
