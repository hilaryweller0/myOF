/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::primitiveMesh

Description
    Additional functions for meshWithCentres to override the cell and face
    centres


\*---------------------------------------------------------------------------*/

#ifndef overrideCentresAreas_H
#define overrideCentresAreas_H

    void overrideFaceCentres(const pointField&);
    void overrideCellCentres(const pointField&);
    void overrideFaceAreas(const vectorField&);
    void overrideCellVols(const scalarList&);
    
    void calcSphFaceCentreAndArea
    (
        const pointField& p,
        vector& fCtr,
        vector& fArea,
        const bool setCentre = true
    ) const;

    void calcSphFaceCentresAndAreas
    (
        const pointField& p, vectorField& fCtrs, vectorField& fAreas,
        const bool setCentres = true
    ) const;
    
    scalar calcSphCellVol
    (
        const vectorField& fCtrs,
        const vectorField& fAreas
    ) const;

    void calcSphCellCentresAndVols
    (
        const vectorField& fCtrs,
        const vectorField& fAreas,
        vectorField& cellCtrs,
        scalarList& cellVols,
        const bool setCentres = true
    ) const;

#endif

// ************************************************************************* //
