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

Application
    setPlaneJet

Description
    Set the initial U and h for a shallow water eqn jet on a beta plane. 
    Also calculate the mesh densities and degrees of freedom

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "cellSet.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cellSet", "cellSetName", "only calculate dofs for a subset of cells"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading initial conditions\n" << endl;

    IOdictionary initDict
    (
        IOobject
        (
            "initialConditions",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Maximum jet velocity
    const dimensionedScalar u0(initDict.lookup("u0"));
    // y location of the centre of the jet
    const dimensionedScalar yc(initDict.lookup("jetCentre"));
    // jet half width
    const dimensionedScalar w(initDict.lookup("jetHalfWidth"));
    // Shallow water height at the equator
    const dimensionedScalar h0(initDict.lookup("h0"));
    
    // Height, location and widths of a bump
    const dimensionedScalar bumpHeight
    (
        initDict.lookupOrDefault<dimensionedScalar>
        (
            "bumpHeight", dimensionedScalar("bumpHeight", dimLength, scalar(0))
        )
    );
    const dimensionedVector bumpCentre
    (
        initDict.lookupOrDefault<dimensionedVector>
        (
            "bumpCentre", dimensionedVector("bumpCentre", dimLength, point::zero)
        )
    );
    const dimensionedVector bumpWidth
    (
        initDict.lookupOrDefault<dimensionedVector>
        (
            "bumpWidth", dimensionedVector("bumpWidth", dimLength, point::zero)
        )
    );
    
    const scalar yb = (yc - w).value();
    const scalar yt = (yc + w).value();

    Info<< "Reading earth properties\n" << endl;

    IOdictionary envDict
    (
        IOobject
        (
            "earthProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // boundary conditions
    wordList hBCs(mesh.boundaryMesh().size(), "empty");
    wordList UBCs(mesh.boundaryMesh().size(), "empty");
    wordList UfBCs(mesh.boundaryMesh().size(), "empty");
    wordList phiBCs(mesh.boundaryMesh().size(), "empty");
    forAll(mesh.boundaryMesh(), patchi)
    {
        const word btype = mesh.boundaryMesh()[patchi].type();
        const word bname = mesh.boundaryMesh()[patchi].name();
        Info << "patch type " << btype << " patch name " << bname << endl;
        if (btype == "wall" || btype == "symmetryPlane")
        {
            hBCs[patchi] = "zeroGradient";
            UBCs[patchi] = "fixedValue";
            UfBCs[patchi] = "fixedValue";
            phiBCs[patchi] = "fixedValue";
        }
        else if (bname == "inlet")
        {
            hBCs[patchi] = "zeroGradient";
            UBCs[patchi] = "calculated";
            UfBCs[patchi] = "calculated";
            phiBCs[patchi] = "fixedValue";
        }
        else if (bname == "outlet")
        {
            hBCs[patchi] = "fixedValue";
            UBCs[patchi] = "calculated";
            UfBCs[patchi] = "calculated";
            phiBCs[patchi] = "calculated";
        }
    }

    // gravity magnitude
    const dimensionedScalar g(envDict.lookup("magg"));
    // beta
    const dimensionedScalar beta(envDict.lookup("beta"));

    dimensionedScalar hmin = h0 - 32/35.*u0*beta/g*w*yc;
    
    volScalarField y("y", mesh.C().component(vector::Y));
    volScalarField yy("yy", (y - yc)/w);
    volScalarField u("u", u0*(1 - 3*sqr(yy) + 3*pow(yy,4) - pow(yy,6)));
    u == max(u, u0*0);

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        vector(1.,0.,0.)*u,
        UBCs
    );
        
    surfaceScalarField yf("yf", mesh.Cf().component(vector::Y));
    surfaceScalarField yyf("yyf", (yf - yc)/w);
    surfaceScalarField uf("uf", u0*(1 - 3*sqr(yyf) + 3*pow(yyf,4) - pow(yyf,6)));
    uf = max(uf, u0*0);

    surfaceVectorField Uf
    (
        IOobject("Uf", runTime.timeName(), mesh),
        vector(1.,0.,0.)*uf,
        UfBCs
    );
    
    // Create initial height and velocity field
    volScalarField htmp
    (
        IOobject("h", runTime.timeName(), mesh),
        mesh,
        h0,
        "zeroGradient"
    );
    
    forAll(htmp, i)
    {
        if (y[i] > yb && y[i] < yt)
        {
            dimensionedScalar hh = h0 - u0*beta/g*w*
            (
                -0.125*w + 16/35.*yc
              + yc*yy[i]*
                (1 - sqr(yy[i]) + 0.6*pow(yy[i],4) - 1/7.*pow(yy[i],6))
              + w*sqr(yy[i])*
                (0.5 - .75*sqr(yy[i]) + 0.5*pow(yy[i],4) - 0.125*pow(yy[i],6))
            );
        
            htmp[i] = hh.value();
        }
        else if(y[i] >= yt)
        {
            htmp[i] = hmin.value();
        }
    }
    htmp.correctBoundaryConditions();
    volScalarField h
    (
        IOobject("h", runTime.timeName(), mesh),
        htmp,
        hBCs
    );
//    volScalarField& h = htmp;
    
    if (bumpHeight.value() > SMALL)
    {
        volScalarField hun("hunperturbed", h);
        hun.write();
        
        volScalarField x("x", mesh.C().component(vector::X));
        
        forAll(h, i)
        {
            if (y[i] > yb && y[i] < yt)
            {
                h[i] += bumpHeight.value()
                       *sqr(Foam::cos(0.5*pi*yy[i]))
                       *Foam::exp(-sqr((x[i]-bumpCentre.value()[0])/bumpWidth.value()[0]))
                       *Foam::exp(-sqr((y[i]-bumpCentre.value()[1])/bumpWidth.value()[1]));
            }
        }
    }
    
    h.write();
    U.write();
    Uf.write();
    volVectorField hU
    (
        IOobject("hU", runTime.timeName(), mesh),
        h*U
    );
    hU.write();
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi", runTime.timeName(), mesh
        ),
        (Uf & mesh.Sf())*fvc::interpolate(h),
        phiBCs
    );
    phi.write();
    
    // List of cells to sum
    labelList sumCells;
    if (args.optionFound("cellSet"))
    {
        const word cellSetName(args.optionRead<string>("cellSet"));
        cellSet cells(mesh, cellSetName);
        sumCells = cells.toc();
    }
    else
    {
        // Select all cells
        sumCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            sumCells[cellI] = cellI;
        }
    }

    Info << "Calculating resolution densities and writing to file " << flush;
    fileName densFile = args.rootPath() / args.caseName() / "meshDensity.dat";
    Info << densFile << endl;
    OFstream os(densFile);
    // min and max x and y of coarse and find resolution regions
    const scalar minxCoarse = 0.5*w.value();
    const scalar maxxCoarse = 1.5*w.value();
    const point fineResCentre(3*w.value(), 1.5*w.value(), 0.5);
    const scalar fineResRadius(0.5*w.value());
    
    label dofCoarse = 0;
    label dofFine = 0;
    scalar areaCoarse = 0;
    scalar areaFine = 0;
    label dofTotal = 0;
    // Loop over the cells
    const pointField& C = mesh.C();
    forAll(C, celli)
    {
        if (C[celli].x() >= minxCoarse && C[celli].x() <= maxxCoarse)
        {
            dofCoarse++;
            areaCoarse += mag(mesh.V()[celli]);
        }
        else if (magSqr(C[celli] - fineResCentre) < sqr(fineResRadius))
        {
            dofFine++;
            areaFine += mag(mesh.V()[celli]);
        }
    }
    areaFine /= mesh.bounds().span().z();
    // loop over the internal faces
    const pointField& Cf = mesh.Cf();
    forAll(Cf, faceI)
    {
        if (Cf[faceI].x() >= minxCoarse && Cf[faceI].x() <= maxxCoarse)
        {
            dofCoarse++;
        }
        else if (magSqr(Cf[faceI] - fineResCentre) < sqr(fineResRadius))
        {
            dofFine++;
        }
    }
    
    // Calculate number of degrees of freedom in the cell set
    forAll(sumCells, i)
    {
        const cell& c = mesh.cells()[sumCells[i]];
        dofTotal += c.size()-2;
    }
    dofTotal *= 0.5;
    dofTotal += sumCells.size();
    
    scalar rhoCoarse = dofCoarse/(3*areaCoarse);
    scalar rhoFine = dofFine/(3*areaFine);
    os << "Coarse region dof = " << dofCoarse << " area = " << areaCoarse
       << " rho = " << rhoCoarse << " dx = " << Foam::sqrt(1./rhoCoarse) << nl;
    os << "Fine region dof = " << dofFine << " area = " << areaFine
       << " rho = " << rhoFine << " dx = " << Foam::sqrt(1./rhoFine) << nl;
    os << "k = " << Foam::sqrt(rhoFine/rhoCoarse) << nl;
    os << "dof = " << dofTotal << nl;

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
