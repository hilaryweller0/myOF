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
    meshAnalysis

Description
    Calculates the non-orthogonality, skewness and some other grid
    characteristics and writes the results relative to geodesic distance
    from a given point

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "TRiSK.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

using namespace Foam;
using namespace constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::validArgs.append("dict");
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const fileName dictName(args.additionalArgs()[0]);
    
    // Read the centre point for the analysis from the dictionary
    IOdictionary dict
    (
        IOobject(dictName, runTime.system(), runTime, IOobject::MUST_READ)
    );
    const point centre(dict.lookup("centre"));
    const point centreHat = unitVector(centre);

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

    volScalarField volRatio
    (
        IOobject("volRatio", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar("volRatio", dimless, scalar(0)),
        "zeroGradient"
    );
    dimensionedScalar Vmean = sum(mesh.V())/mesh.nCells();
    volRatio.internalField() = mesh.V()/Vmean;
    volRatio.correctBoundaryConditions();
    volRatio.write();
    
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

    // Read detHess
    volScalarField detHess
    (
        IOobject("detHess", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT),
        mesh,
        dimensionedScalar("detHess", dimless, scalar(0))
    );

    // Write the mesh analysis for volFields
    fileName ofName
    (
        runTime.rootPath()/runTime.caseName()/runTime.timeName()/"distArea.dat"
    );
    Info << "Writing to " << ofName << endl;
    OFstream ofVol(ofName);
    ofVol << "#dist area nonCentoidaligy detHess\n";
    forAll(mesh.V(), cellI)
    {
        scalar r = mag(mesh.C()[cellI]);
        scalar dist = Foam::acos(point(mesh.C()[cellI]/r) & centreHat);
        scalar area = mesh.V()[cellI]/mesh.depth()[cellI]/sqr(r);
        ofVol << dist << " " << area << " " << nonCentroidalality[cellI]
              << " " << detHess[cellI] << nl;
    }
    
    // Write the mesh analysis for surfaceFields
    fileName ofName2
    (
        runTime.rootPath()/runTime.caseName()/runTime.timeName()/"distMetrics.dat"
    );
    Info << "Writing to " << ofName2 << endl;
    OFstream ofS(ofName2);
    ofS << "#dist dx orthogonality skewness\n";
    for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        scalar r = mag(mesh.Cf()[faceI]);
        scalar dist = Foam::acos(point(mesh.Cf()[faceI]/r) & centreHat);
        
        ofS << dist << " " << dx[faceI]/r << " " << orthogonality[faceI]
            << " " << skewness[faceI] << nl;
    }
    
    // Calculate the scaling factor between the Ringler et al density funciton
    // and the cell areas
    if (dict.found("alpha"))
    {
        scalar alpha = readScalar(dict.lookup("alpha"));
        scalar beta  = readScalar(dict.lookup("beta"));
        scalar gamma = readScalar(dict.lookup("gamma"));
    
        scalar sum = 0;
        forAll(mesh.C(), cellI)
        {
            scalar dist = Foam::acos(point(unitVector(mesh.C()[cellI])) & centreHat);
            scalar rho = 0.5/(1+gamma)*(Foam::tanh((beta-dist)/alpha)+1) + gamma;
            sum += 1/Foam::sqrt(rho);
        }
        
        const scalar theoryScale = 4*pi/sum;
        Info << "Scaling factor between theoretical cell area and achieved cell "
             << "area = " << theoryScale << endl;
        
        // Output the theoretical dx and cell area
        OFstream ofVol
        (
            runTime.rootPath()/runTime.caseName()
           /runTime.timeName()/"distAreaTheory.dat"
        );
        ofVol << "#dist dx area" << endl;
        forAll(mesh.C(), cellI)
        {
            scalar r = mag(mesh.C()[cellI]);
            scalar dist = Foam::acos(point(mesh.C()[cellI]/r) & centreHat);
            scalar rho = 0.5/(1+gamma)*(Foam::tanh((beta-dist)/alpha)+1) + gamma;
            scalar area = theoryScale/Foam::sqrt(rho);
            // Cell-centre to cell centre distance for a regular hexagon
            scalar dx = Foam::sqrt(area*2/3*Foam::tan(pi/3));
            ofVol << dist << " " << dx << " " << area << nl;
        }
    }
}


// ************************************************************************* //

