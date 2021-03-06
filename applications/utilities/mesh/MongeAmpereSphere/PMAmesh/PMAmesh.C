/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    PMAmesh

Description
    Solves the Parabolic Monge-Ampere equations to move a mesh. 
    Reads in either a list of points for the initial cell centres or cell
    centres from an existing OpenFOAM case.

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"
#include "fvCFD.H"
#include "InitialPoints.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName()
         <<  " region " << meshRegion << endl;

    // Create the mesh of the sphere based on the regoin "meshRegion"
    fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        ),
        false,
        fvMeshWithDual::SPHERICALDIST
    );
    
    // Open control dictionary
    IOdictionary controlDict
    (
        IOobject
        (
            args.executable() + "Dict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    Info << "Reading initial points and required resolution" << endl;
    autoPtr<InitialPoints> initialPoints(InitialPoints::New(controlDict));
    const dimensionedScalar& earthRadius = mesh.earthRadius();

    // Create the mesh displacement potential, Phi
    volScalarField Phi
    (
        IOobject("Phi",runTime.timeName(),mesh,IOobject::NO_READ),
        mesh,
        dimensionedScalar("Phi", dimArea, scalar(0))
    );
    Phi.oldTime();
    
    volVectorField rHat = mesh.C()/mag(mesh.C());
    
    volVectorField gradPhi("gradPhi", fvc::grad(Phi));
    
    volScalarField detHess
    (
        IOobject("detHess", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(1))
    );
    volTensorField Hessian("Hessian", fvc::grad(gradPhi));
    volScalarField del2Phi("del2Phi", fvc::laplacian(Phi));
    
    volScalarField monitor
    (
        IOobject("monitor", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("", dimless, scalar(0))
    );
    forAll(monitor, cellI)
    {
        monitor[cellI] = 1./
        (
            initialPoints().requiredResolution(mesh.C()[cellI])
        );
    }
    volScalarField equidist("equidist", monitor*detHess);
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << " " << flush;
        
        dimensionedScalar equiDistMean
           = sum((detHess - del2Phi)*mesh.V())
            /sum(1./monitor*mesh.V());
            
        fvScalarMatrix PhiEqn
        (
            fvm::laplacian(Phi)
          - equiDistMean/monitor
          - del2Phi
          + detHess
        );
        PhiEqn.setReference(0, scalar(0));
        PhiEqn.solve();

        gradPhi = fvc::grad(Phi);
        Hessian = fvc::grad(gradPhi);
        Hessian = 0.5*(Hessian + Hessian.T());
        del2Phi = fvc::laplacian(Phi);
        forAll(detHess, cellI)
        {
            detHess[cellI] = det
            (
                diagTensor::one + Hessian[cellI]
            );
        }
        // Correct the linear terms of detHess
        detHess += del2Phi - tr(Hessian);
        
        Info << "Time = " << runTime.timeName()
             << " determinant goes from " << min(detHess).value()
             << " to " << max(detHess).value()
             << " equidist goes from " << min(equidist).value()
             << " to " << max(equidist).value() << endl;
        forAll(monitor, cellI)
        {
            monitor[cellI] = 1./
            (
                initialPoints().requiredResolution
                (
                    earthRadius.value()*unitVector(mesh.C()[cellI] + gradPhi[cellI])
                )
            );
        }
        
        equidist = monitor*detHess;
        
        if (min(detHess).value() < 0)
        {
            runTime.writeAndEnd();
        }
    }
    
    Phi.write();
    gradPhi.write();
    monitor.write();
    equidist.write();
    detHess.write();
    Hessian.write();
    del2Phi.write();
    
    pointField newPoints
    (
        earthRadius.value()*unitVector(mesh.C() + gradPhi)
    );
    OFstream os("newPoints");
    os << newPoints;
    //newPoints.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
