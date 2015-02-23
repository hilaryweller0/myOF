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
#include "IFstream.H"
#include "OFstream.H"
#include "dimensionedTypes.H"
//#include "polarPoint.H"
#include "Time.H"
//#include "IOmanip.H"
//#include "findInterpWeights.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * ** Member functions  * * * * * * * * * * * * *  //



//Foam::scalar Foam::VoronoiSphereMeshing::controls::interp
//(
//    const std::map<Point, K::FT, K::Less_xyz_3>& f,
//    const Point& p
//) const
//{
//    typedef std::vector<std::pair<Point, K::FT> > Point_coordinate_vector;
//    typedef CGAL::Data_access<std::map<Point, K::FT, K::Less_xyz_3> > Value_access;

//    Vector_3 normal(p-CGAL::ORIGIN);
//    
//    Point_coordinate_vector coords;
//    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool>
//        result = CGAL::surface_neighbor_coordinates_3
//    (
//        origMeshTri,
//        p,
//        normal,
//        std::back_inserter(coords)
//    );
//    
//    if (!result.third)
//    {
//        FatalErrorIn("Foam::VoronoiSphereMeshing::controls::interp")
//            << " the coordinate computation was not successful onto point ("
//            << p.x() << " " << p.y() << " " << p.z() << ')' << exit(FatalError);
//    }
//    
//    K::FT norm = result.second;

//    // Linear interpolation
//    return CGAL::linear_interpolation
//    (
//        coords.begin(),
//        coords.end(),
//        norm,
//        Value_access(f)
//    );
//}

// ************************************************************************* //
