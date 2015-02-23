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

Typedefs
    CGALTriangulation3Ddefs

Description
    CGAL data structures used for 3D Delaunay meshing.

    Define CGAL_INEXACT to use Exact_predicates_inexact_constructions kernel
    otherwise the more robust but much less efficient
    Exact_predicates_exact_constructions will be used.

    Define CGAL_HIERARCHY to use hierarchical Delaunay triangulation which is
    faster but uses more memory than the standard Delaunay triangulation.

\*---------------------------------------------------------------------------*/

#ifndef CGALTriangulation3Ddefs_H
#define CGALTriangulation3Ddefs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CGAL/Delaunay_triangulation_3.h"
#include "CGAL/Triangulation_cell_base_with_circumcenter_3.h"

#include "indexedVertex.H"
#include "indexedCell.H"

#ifdef CGAL_INEXACT
    // Fast kernel using a double as the storage type but the triangulation
    // may fail
    #include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#else
    // Very robust but expensive kernel
    #include "CGAL/Exact_predicates_exact_constructions_kernel.h"
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#endif

typedef CGAL::indexedVertex<K>                           Vb;
typedef CGAL::indexedCell<K>                             Cb;

#ifdef CGAL_HIERARCHY
    // Data structures for hierarchical Delaunay triangulation which is more
    // efficient but also uses more storage
    #include "CGAL/Triangulation_hierarchy_3.h"
    typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
    typedef CGAL::Triangulation_data_structure_3<Vbh, Cb>    Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds>           Triangulation;
    typedef CGAL::Triangulation_hierarchy_3<Triangulation>   HTriangulation;
#else
    // Data structures for standard Delaunay triangulation
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb>     Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds>           Triangulation;
    typedef Triangulation                                    HTriangulation;
#endif

#include <CGAL/surface_neighbor_coordinates_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/spatial_sort.h>

typedef K::Vector_3 Vector_3;
typedef CGAL::Tetrahedron_3<K> Tet;
typedef CGAL::Triangle_3<K> Triangle;
typedef CGAL::Line_3<K> Line;

//#include <CGAL/Surface_mesh_default_triangulation_3.h>
//#include <CGAL/Complex_2_in_triangulation_3.h>
//#include <CGAL/make_surface_mesh.h>
//#include <CGAL/Implicit_surface_3.h>

//// default triangulation for Surface_mesher
//typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

//// c2t3
//typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

//typedef Tr::Geom_traits GT;
//typedef GT::Sphere_3 Sphere_3;
//typedef GT::Point_3 Point_3;
//typedef GT::FT FT;

//typedef FT (*Function)(Point_3);

//typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

//inline FT sphereFunction (Point_3 p)
//{
//    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
//    return x2+y2+z2-1;
//}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
