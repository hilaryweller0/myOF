/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Application
    orthogonality

Description
    Calculates the non-orthogonality, skewness and some other grid
    characteristics between the primal and dual grids and write
    surfaceScalarFields on the primal and dual with the results

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "TRiSK.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createMeshWithDual.H"
    Foam::fvMeshWithDual centroidalMesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        false,
        fvMeshWithDual::SPHERICALDIST
    );
    #include "createDualMesh.H"
    centroidalMesh.setDual(dualMesh);
    
    Foam::fvMeshWithDual centroidalDualMesh
    (
        Foam::IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        false,
        fvMeshWithDual::SPHERICALDIST
    );
    
    surfaceScalarField orthogonality
    (
        IOobject("orthogonality", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("degrees", dimless, scalar(0))
    );
    
    surfaceScalarField skewness
    (
        IOobject("skewness", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("degrees", dimless, scalar(0))
    );
    
    forAll(mesh.Cf(), faceI)
    {
        const point& Co = mesh.C()[mesh.owner()[faceI]];
        const point& Cn = mesh.C()[mesh.neighbour()[faceI]];
        label faced = mesh.dualFaceMap()[faceI];
        const point& Cdo = dualMesh.C()[dualMesh.owner()[faced]];
        const point& Cdn = dualMesh.C()[dualMesh.neighbour()[faced]];
    
        vector deltaP = Co - Cn;
        vector deltaD = Cdo - Cdn;

        orthogonality[faceI] = Foam::acos
        (
            (deltaP & deltaD)/(mag(deltaP)*mag(deltaD))
        )*180./pi;
        orthogonality[faceI] = mag(90-orthogonality[faceI]);
        
        skewness[faceI] = angleBetween(mesh.Cf()[faceI], 0.5*(Cdo + Cdn))
                         /angleBetween(Co, Cn);
    }
    
    surfaceScalarField dx("dx", 1./mesh.deltaCoeffs());
    dx.write();
    
    orthogonality.write();
    skewness.write(); 
    
    // Dual versions
    surfaceScalarField orthogD
    (
        IOobject("orthogonality", dualMesh.time().timeName(), dualMesh),
        mesh.dualMap(orthogonality)
    );
    orthogD.write();
    surfaceScalarField skewD
    (
        IOobject("skewness", dualMesh.time().timeName(), dualMesh),
        mesh.dualMap(skewness)
    );
    skewD.write();
    
    const dimensionedScalar Vtot = sum(mesh.V());

    Info << "n cells = " << mesh.nCells()
         << " ndofs = " << mesh.nCells()+mesh.nInternalFaces()
         << "\nmean dx = " << (TRiSK::domainIntegrate(dx)/Vtot).value()
         << " max/min dx = " << (max(dx)/min(dx)).value()
         << " max dx = " << (max(dx)).value()
         << "\nmax non-orthogonality = " << max(orthogonality).value()
         << " mean non-orthogonality = "
         << (TRiSK::domainIntegrate(orthogonality)/Vtot).value()
         << "\nmax skewness = " << max(skewness).value()
         << " mean skewness = "
         << (TRiSK::domainIntegrate(skewness)/Vtot).value() << endl;
    
    // For calculating the distance between the cell centroid and the cell centre
    // we must create a new mesh with the cell centres at the centroids
    
    volScalarField nonCentroidalality
    (
        IOobject("nonCentroidalality", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("nonCentroidalality", dimless, scalar(0))
    );
    
    nonCentroidalality.internalField() = mag(mesh.C()[0])*
                       mag(unitVector(mesh.C()) - unitVector(centroidalMesh.C()))
                     / Foam::sqrt(mesh.V()/mesh.depth().internalField());

    Info << "max nonCentroidalality = " << max(nonCentroidalality).value()
         << " mean nonCentroidalality = "
         << (fvc::domainIntegrate(nonCentroidalality)/Vtot).value() << endl;
    
    volScalarField nonCentroidalityDual
    (
        IOobject("nonCentroidality", dualMesh.time().timeName(), dualMesh),
        dualMesh,
        dimensionedScalar("nonCentroidality", dimless, scalar(0))
    );
    
    nonCentroidalityDual.internalField() = mag(dualMesh.C()[0])*
                       mag(unitVector(dualMesh.C()) - unitVector(centroidalDualMesh.C()))
                     / Foam::sqrt(dualMesh.V()/dualMesh.depth().internalField());

    Info << "Dual max nonCentroidalality = " << max(nonCentroidalityDual).value()
         << " mean nonCentroidalality = "
         << (fvc::domainIntegrate(nonCentroidalityDual)/Vtot).value() << endl;
}


// ************************************************************************* //

