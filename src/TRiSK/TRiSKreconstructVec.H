/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

InNamespace
    Foam::TRiSK

Description
    Reconstruct the vector field at faces from the normal components

SourceFiles
    TRiSKreconstructVec.C

\*---------------------------------------------------------------------------*/


#ifndef TRiSKreconstructVec_H
#define TRiSKreconstructVec_H

#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace TRiSK functions Declaration
\*---------------------------------------------------------------------------*/

namespace TRiSK
{
    //- TRiSK reconstruction from primal fluxes to maintain divergence on the dual
    tmp<surfaceVectorField> divReconstruct(const surfaceScalarField& phi);

    //- TRiSK reconstruction to maintain divergence on the dual
    tmp<surfaceVectorField> divReconstruct(const tmp<surfaceScalarField>&);

    //- TRiSK reconstruction from primal fluxes to maintain divergence on the dual
    tmp<surfaceVectorField> divReconstruct3d(const surfaceScalarField& phi);

    //- TRiSK reconstruction to maintain divergence on the dual
    tmp<surfaceVectorField> divReconstruct3d(const tmp<surfaceScalarField>&);

    //- Reconstruct cell centre velocity from face velocites in direction d
    tmp<volVectorField> reconstruct
    (
        const surfaceScalarField& v,
        const surfaceVectorField& d
    );
    
    //- Reconstruct face veloctiy from face velocities in direction d
    tmp<surfaceVectorField> faceReconstruct
    (
        const surfaceScalarField& v,
        const surfaceVectorField& d
    );
    
    //- map velocity component from dual in dirction j to primal in direction i
    tmp<surfaceScalarField> speedMap
    (
        const surfaceScalarField& v,
        const surfaceVectorField& j,
        const surfaceVectorField& i
    );
    
    //- map fluxes from dual in directio j to primal in direction i
    tmp<surfaceScalarField> fluxMap
    (
        const surfaceScalarField& vS,
        const surfaceVectorField& j,
        const surfaceVectorField& i
    );
    
    //- map the circulation on the dual to the primal flux
    tmp<surfaceScalarField> circToFlux(const surfaceScalarField& vS);
    
    // The off-diaganol of the H operator (mapping fluxes from direction j
    // on the dual to direction i on the primal
    tmp<surfaceScalarField> circToFluxOffDiag(const surfaceScalarField& vS);
    
    //- map the velocity in the d direction to the normal direction
    tmp<surfaceScalarField> ddirToFlux(const surfaceScalarField& v);
    
    // The off-diaganol of the H operator (mapping velocities from
    // direction d on the primal to direction i on the primal
    tmp<surfaceScalarField> ddirToFluxOffDiag(const surfaceScalarField& v);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
