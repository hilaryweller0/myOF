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

Foam::VoronoiSphereMeshing::Edge Foam::VoronoiSphereMeshing::edgeOpposite
(
    const Cell_handle& ch,
    const Vertex_handle& v0,
    const Vertex_handle& v1
) const
{
    // Find the indices of the far points
    int i,j;
    Cell c = *ch;
    c.has_vertex(v0, i);
    c.has_vertex(v1,j);
    int k = 0;
    if (k == min(i,j))
    {
        if (max(i,j) == 1) k = 2;
        else k = 1;
    }
    int l = 6 - i - j - k;

    return Edge(ch, k, l);
}
        

Foam::VoronoiSphereMeshing::Edge Foam::VoronoiSphereMeshing::edgeBetween
(
    const Vertex_handle& v0, const Vertex_handle& v1
) const
{
    Cell_handle c;
    int i, j;
    if (!is_edge(v0, v1, c, i, j))
    {
        FatalErrorIn("VoronoiSphereMeshing::edgeBetween")
            << " no edge in between vertices at points " << topoint(v0->point())
            << " and " << topoint(v1->point())
            << exit(FatalError);
    }
    return Edge(c,i,j);
}

Foam::point Foam::VoronoiSphereMeshing::circumcenter
(
    const Cell_handle& c
) const
{
    const Vertex_handle& v = vCentral();

    Vertex_handle v0 = c->vertex(0);
    Vertex_handle v1 = c->vertex(1);
    Vertex_handle v2 = c->vertex(2);

    if      (v0 == v) v0 = c->vertex(3);
    else if (v1 == v) v1 = c->vertex(3);
    else if (v2 == v) v2 = c->vertex(3);
    else if (v != c->vertex(3))
    {
        FatalErrorIn("VoronoiSphereMeshing::circumcenter") << " vertex at point "
            << topoint(v->point()) << " is not a vertex of the cell"
            << exit(FatalError);
    }

    return topoint(construct_circumcenter(v0->point(), v1->point(), v2->point()));
}

bool Foam::VoronoiSphereMeshing::sphericalCell(const Cell_handle& c) const
{
    const Vertex_handle& v = vCentral();
    
    if
    (
        c->vertex(0) != v && c->vertex(1) != v
     && c->vertex(2) != v && c->vertex(3) != v
    )
    {
        return false;
    }
    return true;
}

Foam::Pair<Foam::point> Foam::VoronoiSphereMeshing::pointsOpposite(const Edge& e) const
{
    Foam::Pair<Foam::point> pts;
    Foam::Pair<bool> found(false, false);
    
    // The two vertices connected by this edge
    Vertex_handle vi = e.first->vertex(e.second);
    Vertex_handle vj = e.first->vertex(e.third);

    Cell_circulator ced = incident_cells(e);
    Cell_circulator cedStart = ced;
    do
    {
        if (!found[0] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[0]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    pts[0] = topoint(v->point());
                    found[0] = true;
                }
            }
        }
        else if (!found[1] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[1]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    pts[1] = topoint(v->point());
                    found[1] = true;
                }
            }
        }
    } while (++ced != cedStart);
    
    return pts;
}

Foam::Pair<Foam::label> Foam::VoronoiSphereMeshing::verticesOpposite(const Edge& e) const
{
    Foam::Pair<label> vs;
    Foam::Pair<bool> found(false, false);
    
    // The two vertices connected by this edge
    Vertex_handle vi = e.first->vertex(e.second);
    Vertex_handle vj = e.first->vertex(e.third);

    Cell_circulator ced = incident_cells(e);
    Cell_circulator cedStart = ced;
    do
    {
        if (!found[0] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[0]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    vs[0] = v->index();
                    found[0] = true;
                }
            }
        }
        else if (!found[1] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[1]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    vs[1] = v->index();
                    found[1] = true;
                }
            }
        }
    } while (++ced != cedStart);
    
    return vs;
}

void Foam::VoronoiSphereMeshing::neighbourCells
(
    Foam::Pair<Cell_handle>& neibCells,
    Foam::Pair<Foam::point>& pOut,
    const Edge& e
) const
{
    Foam::Pair<bool> found(false, false);
    
    // The two vertices connected by this edge
    Vertex_handle vi = e.first->vertex(e.second);
    Vertex_handle vj = e.first->vertex(e.third);

    Cell_circulator ced = incident_cells(e);
    Cell_circulator cedStart = ced;
    do
    {
        if (!found[0] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[0]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    pOut[0] = topoint(v->point());
                    neibCells[0] = ced;
                    found[0] = true;
                }
            }
        }
        else if (!found[1] && !is_infinite(ced))
        {
            for(label iv = 0; iv < 4 && !found[1]; iv++)
            {
                Vertex_handle v = ced->vertex(iv);
                if (v != vCentral() && !is_infinite(v) && v != vi && v != vj)
                {
                    pOut[1] = topoint(v->point());
                    neibCells[1] = ced;
                    found[1] = true;
                }
            }
        }
    } while (++ced != cedStart && !found[1]);
    
    if (!found[1])
    {
        FatalErrorIn("VoronoiSphereMeshing::neighbourCells")
            << " cannot find 2 finite cells incident to edge with points at\n"
            << topoint(e.first->vertex(e.second)->point()) << " and\n"
            << topoint(e.first->vertex(e.third)->point())
            << exit(FatalError);
    }
}

void Foam::VoronoiSphereMeshing::incidentVertices
(
    const Vertex_handle& vit,
    std::list<Vertex_handle>& incidentVerts
) const
{
    // Ciculate around the cells around the vertex and select the far points
    Cell_circulator ccirc = incident_cells(edgeBetween(vit,vCentral()));
    const Cell_circulator cStart = ccirc;
    Edge e = edgeOpposite(ccirc, vit, vCentral());
    Vertex_handle v0 = ccirc->vertex(e.second);
    Vertex_handle v1 = ccirc->vertex(e.third);
    
    for (ccirc++; ccirc != cStart; ccirc++)
    {
        e = edgeOpposite(ccirc, vit, vCentral());
        Vertex_handle v3 = ccirc->vertex(e.second);
        Vertex_handle v4 = ccirc->vertex(e.third);
    
        if (v3 == v0 || v3 == v1)
        {
            incidentVerts.push_back(v3);
            v0 = v3;
            v1 = v4;
        }
        else
        {
            incidentVerts.push_back(v4);
            v0 = v4;
            v1 = v3;
        }
    }
    incidentVerts.push_back(v1);
}


// ************************************************************************* //
