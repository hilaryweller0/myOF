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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef CGAL_INEXACT

inline Foam::VoronoiSphereMeshing::pointFromPoint Foam::VoronoiSphereMeshing::topoint
(
    const Point& P
) const
{
    return reinterpret_cast<pointFromPoint>(P);
}

inline Foam::VoronoiSphereMeshing::PointFrompoint Foam::VoronoiSphereMeshing::toPoint
(
    const Foam::point& p
) const
{
    return reinterpret_cast<PointFrompoint>(p);
}

#else

inline Foam::VoronoiSphereMeshing::pointFromPoint Foam::VoronoiSphereMeshing::topoint
(
    const Point& P
) const
{
    return Foam::point
    (
        CGAL::to_double(P.x()),
        CGAL::to_double(P.y()),
        CGAL::to_double(P.z())
    );
}

inline Foam::VoronoiSphereMeshing::PointFrompoint Foam::VoronoiSphereMeshing::toPoint
(
    const Foam::point& p
) const
{
    return Point(p.x(), p.y(), p.z());
}

#endif

inline Foam::VoronoiSphereMeshing::PointFrompoint Foam::VoronoiSphereMeshing::toPointWithSize
(
    const Foam::point& p, const scalar r
) const
{
    return toPoint(p*r/mag(p));
}


// ************************************************************************* //
