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

#include "VoronoiSphereMeshing.H"

// * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * //

Foam::label boolSum(const Foam::boolList& a)
{
    Foam::label s = 0;
    forAll(a, i)
    {
        s += a[i];
    }
    return s;
}

// * * * * * * * * * * * * * ** Member functions  * * * * * * * * * * * * *  //

Foam::vector Foam::VoronoiSphereMeshing::springForce
(
    const Foam::point& x,
    const Foam::point& y,
    const label np
) const
{
    Foam::point midPoint = 0.5*(x + y);
    vector delta = y - x;
    scalar magDelta = mag(delta);
    delta *= (np*magDelta - initialPoints.requiredResolution(midPoint))/magDelta;
    return delta;
}

void Foam::VoronoiSphereMeshing::whichPointsOk()
{
    meshOk = boolList(number_of_vertices(), !smoothAll);

    if (!smoothAll)
    {
        // vertices in a non-uniform region of mesh need smoothing
        label iv = 1;
        Finite_vertices_iterator vih = finite_vertices_begin();
        vih++;
        do
        {
            if (meshOk[iv])
            {
                // Current point
                Foam::point defVert0 = topoint(vih->point());
                // The min and max resolution at this vertex
                scalar dx2Min = -1;
                scalar dx2Max = -1;

                // All vertices connected to current vertex
                std::list<Vertex_handle> incidentVerts;
                incident_vertices(vih, std::back_inserter(incidentVerts));
                // loop around the incident vertices to calculate distances
                for
                (
                    std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
                    vin != incidentVerts.end();
                    ++vin
                )
                {
                    // only for vertices on the sphere
                    if (*vin != infinite_vertex() && *vin != vCentral())
                    {
                        Foam::point otherPoint = topoint((*vin)->point());
                        scalar dx2 = magSqr(defVert0 - otherPoint);

                        if (dx2Min < 0)
                        {
                            dx2Min = dx2;
                            dx2Max = dx2;
                        }
                        else
                        {
                            if (dx2 < dx2Min) dx2Min = dx2;
                            else if (dx2 > dx2Max) dx2Max = dx2;
                        }
                    }
                }
                // smooth vertex and its neighbours if large aspect ratio
                // or not hexagonal
                if (dx2Max/dx2Min > 1.5) // || incidentVerts.size() != 8)
                {
                    meshOk[iv] = false;

                }
            }
            iv++;
        } while(++vih != finite_vertices_end());

        Info << label(number_of_vertices()-boolSum(meshOk)) << " of "
            << label(number_of_vertices())-1 << " vertices will be smoothed" << endl;
    }
}


void Foam::VoronoiSphereMeshing::reIndex()
{
    label i = 1;
    Finite_vertices_iterator vit = finite_vertices_begin();
    vit++;
    do
    {
        vit->index() = i++;
    } while(++vit != finite_vertices_end());
}


void Foam::VoronoiSphereMeshing::movePoint
(
    const Vertex_handle& vh,
    const Point& p
)
{
    // Remember an incident vertex to start and remember the vertex index
    Cell_handle ch = vh->cell();
    Vertex_handle old_neighbor = ch->vertex(ch->index(vh) == 0 ? 1 : 0);
    label nv = vh->index();

    //move_point(vh, p);
    //vh->index() = nv;

    remove(vh);
    insert(p, old_neighbor->cell())->index()=nv;
}

Foam::point Foam::VoronoiSphereMeshing::intersection
(
    const Foam::point& a1, const Foam::point& b1
) const
{
    const Foam::point a = 0.9*a1;
    const Foam::point b = 0.9*b1;

    // point to be found
    Foam::point intersect;

    // First locate the cell containing a
    const Cell_handle teta = locate(toPoint(a));
    
    // find the triangular facet of teta of which b is the other side
    bool facetFound = false;
    for(int i = 0; i < 4 && !facetFound; i++)
    {
        const Cell_handle tetb = teta->neighbor(i);
        Locate_type loc;
        int j, k;
        if (side_of_cell(toPoint(b), tetb, loc, j, k) == 1)
        {
            facetFound = true;
            const Point t0 = teta->vertex(i > 0 ? 0 : 1)->point();
            const Point t1 = teta->vertex(i > 1 ? 1 : 2)->point();
            const Point t2 = teta->vertex(i > 2 ? 2 : 3)->point();
            
            Object interObj = CGAL::intersection
            (
                Triangle(t0,t1,t2),
                Line(toPoint(a), toPoint(b))
            );

            if (const Point * p = CGAL::object_cast<Point>(&interObj))
            {
                intersect = topoint(*p);
            }
            else
            {
                FatalErrorIn("VoronoiSphereMeshing::intersection")
          << " cannot find intersection between the line between the points\n"
                    << a << " and\n" << b << nl
                    << "and the triangle with points\n"
                    << topoint(t0) << nl
                    << topoint(t1) << nl
                    << topoint(t2) << exit(FatalError);
            }
        }
    }
    if (!facetFound)
    {
        FatalErrorIn("VoronoiSphereMeshing::intersection")
            << " cannot find a facet of teta which intersects points\n"
            << a << " and\n" << b << exit(FatalError);
    }
    
    return unitVector(intersect)*mag(a1);
}


void Foam::VoronoiSphereMeshing::insertPoints(const pointField& points)
{
    label nVert = number_of_vertices();

    if (nVert == 0)
    {
        insert(Point(0.,0.,0.))->index() = nVert++;
    }

    Vertex_handle vh = finite_vertices_begin();

    // Add the points and index them
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        insert(toPoint(unitVector(p)*radius), vh->cell())->index() = nVert++;
        vh++;
    }
}

void Foam::VoronoiSphereMeshing::insertPoints
(
    const pointField& points,
    List<Vertex_handle>& vh
)
{
    if (points.size() > vh.size())
    {
        FatalErrorIn("VoronoiSphereMeshing::insertPoints(const pointField& points, List<Vertex_handle>& vh") << " there are more points to insert than vertices to insert them near but points.size() = " << points.size() << " and vh.size() = " << vh.size() << exit(FatalError);
    }

    label nVert = number_of_vertices();

    if (nVert == 0)
    {
        insert(Point(0.,0.,0.))->index() = nVert++;
    }

    // Add the points and index them
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        insert(toPoint(unitVector(p)*radius),vh[i]->cell())->index()= nVert++;
    }
}


Foam::label Foam::VoronoiSphereMeshing::removePoints(const boolList& removePoint)
{
    label nPointsRemove = 0;
    
    std::list<Vertex_handle>  vertices;
    incident_vertices(infinite_vertex(), std::back_inserter(vertices));
    
    for
    (
        std::list<Vertex_handle>::iterator v_set_it = vertices.begin();
        v_set_it != vertices.end();
        v_set_it++
    )
    {
        if (removePoint[(*v_set_it)->index()])
        {
            remove(*v_set_it);
            nPointsRemove++;
        }
    }

//    Finite_vertices_iterator vit = finite_vertices_begin();
//    vit++;
//    do
//    {
//        if (removePoint[vit->index()] == true)
//        {
//            Finite_vertices_iterator vir = vit;
//            vit--;
//            remove(vir);
//            nPointsRemove++;
//        }
//    } while (++vit != finite_vertices_end());

    if (nPointsRemove > 0 && number_of_vertices() > 4) reIndex();
    
    if (number_of_vertices() > 1 && number_of_vertices() <= 4)
    {
        WarningIn("VoronoiSphereMeshing::removePoints")
            << "only " << label(number_of_vertices())
            << " points left" << endl;
    }

    return nPointsRemove;
}


// ************************************************************************* //
