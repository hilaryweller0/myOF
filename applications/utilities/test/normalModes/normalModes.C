/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

Application
    normalModes

Description
    Calculates the normal modes, amplitudes and frequencies of a shallow
    water model on the sphere

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "Time.H"
#include "meshWithDual.H"
#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"
#include "OFstream.H"

#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace Foam;
//using namespace std;
//using namespace Eigen;

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("shallowWaterCode");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshWithDual.H"
    #include "createFields.H"
    fileName shallowWaterCode = args.additionalArgs()[0];

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // The model response matrix
    Eigen::MatrixXd M(h.size() + Uf.size(), h.size() + Uf.size());

    forAll(h, cellI)
    {
        runTime.setTime(runTime.startTime(), runTime.startTimeIndex());
        Info << "Creating model response matrix for cell: " << cellI
             << " out of " << mesh.nCells() << endl;
        h.internalField() = 0;
        Uf.internalField() = vector(0,0,0);
        h[cellI] = 1;

        h.write();
        Uf.write();
        
        bool err = system(shallowWaterCode.c_str());
        if (err) FatalErrorIn(shallowWaterCode) << exit(FatalError);
        
        Info << "Reading h and Uf for model response matrix" << endl;
        runTime++;
        h = volScalarField
        (
            IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        Uf = surfaceVectorField
        (
            IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        
        forAll(h, i) { M(cellI, i) = h[i]; }
        forAll(Uf, i) { M(cellI, h.size()+i) = (Uf[i] & mesh.ddir()[i]); }
    }
    forAll(Uf, faceI)
    {
        runTime.setTime(runTime.startTime(), runTime.startTimeIndex());
        Info << "Creating model response matrix for face: " << faceI
             << " out of " << mesh.nInternalFaces() << endl;
        h.internalField() = 0;
        Uf.internalField() = vector(0,0,0);
        Uf[faceI] = mesh.ddir()[faceI];

        h.write();
        Uf.write();

        bool err = system(shallowWaterCode.c_str());
        if (err) FatalErrorIn(shallowWaterCode) << exit(FatalError);

        Info << "Reading h and Uf for model response matrix" << endl;
        runTime++;
        h = volScalarField
        (
            IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        Uf = surfaceVectorField
        (
            IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );

        forAll(h, i) { M(h.size() + faceI, i) = h[i]; }
        forAll(Uf, i)
        {
            M(h.size() + faceI, h.size()+i) = (Uf[i] & mesh.ddir()[i]);
        }
    }
    
    Info << "\nChecking if matrix is self-adjoint" << endl;
    bool selfAdj = M.isApprox(M.adjoint());
    Info << "M selfAdj = " << selfAdj << endl;
    
    Info << "\nCalculating eigenvalues" << endl;
    
    Eigen::EigenSolver<MatrixXd> eigensolver(M);
    Eigen::VectorXcd eigVals = eigensolver.eigenvalues();
    Eigen::VectorXcd eigValsTest = eigVals;
    Eigen::MatrixXcd eigVecs = eigensolver.eigenvectors();

    scalarField amplification(M.rows());
    scalarField freq(M.rows());
    
    OFstream evs("eigenValues");
    forAll(freq, i)
    {
        amplification[i] = abs(eigVals[i]);
        freq[i] = arg(eigVals[i])/runTime.deltaT().value();
        Info << "Eigen value " << std::real(eigVals[i]) << ", "
             << std::imag(eigVals[i]) << endl;
        evs << std::real(eigVals[i]) << " " << std::imag(eigVals[i]) << nl;
    }
    
    {
        OFstream ff("freq");
        ff << freq;
        OFstream af("amplification");
        af << amplification;
        
    }
    
    surfaceScalarField v = Uf & mesh.ddir();
    Info<< "\nStarting time loop for writing out eigenvectors\n" << endl;
    runTime.setTime(runTime.startTime(), runTime.startTimeIndex());    
    Info << "\nTime = ";
    for (label ivec = 0; ivec < M.rows(); runTime+= scalar(1), ivec++)
    {
        Info << runTime.timeName() << " ";
        forAll(h, cellI) h[cellI]  = std::real(eigVecs.col(ivec)[cellI]);
        forAll(v,faceI) v[faceI] = std::real(eigVecs.col(ivec)[h.size()+faceI]);
        Uf = Hops::faceReconstruct(v*mesh.magSf(), mesh.ddir());
        h.write();
        Uf.write();
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
