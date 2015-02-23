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
    MongeAmpereSphere

Description
    Solves the Monge-Ampere equations to move a mesh on the surface of a
    sphere. Monitor function defined analytically.

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
#include "VectorSpaceFunctions.H"
#include "monitorFunction.H"
#include "meshToEdges.H"
#include "polyFit.H"
#include "faceToPointReconstruct.H"
#include "faceToPoint.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    // Check for plotting non-default region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName()
         <<  " region " << meshRegion << endl;

    // Create the mesh of the sphere based on the regoin "meshRegion"
    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    
    // Create the mesh to be moved
    fvMesh rMesh
    (
        Foam::IOobject
        (
            "rMesh", runTime.timeName(), runTime,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        )
    );

    // maps for mapping from cell centres to internal edges
    meshToEdges<polyFit<1> > mte(mesh);
    //faceToPoint<polyFit<1> > ftp(mesh);
    
    // Open control dictionary
    IOdictionary controlDict
    (
        IOobject
        (
            args.executable() + "Dict", runTime.system(), runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    // The monitor funciton
    autoPtr<monitorFunction> monitorFunc(monitorFunction::New(controlDict));
    
    // Read the ratio of the Hessian to volume to use
    const scalar HessianVolumeRatio
    (
        readScalar(controlDict.lookup("HessianVolumeRatio"))
    );

    dimensionedScalar Vtot = sum(mesh.V());

    // boost the Laplacian to stabilise
    scalar boostLaplacian(readScalar(controlDict.lookup("boostLaplacian")));

    // Use a wide or standard stencil
    const Switch wideStencil(controlDict.lookup("wideStencil"));

    #include "createFields.H"

    // Use time-steps instead of iterations to solve the Monge-Ampere eqn
    bool converged = false;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << flush << nl;

        boostLaplacian = max
        (
            boostLaplacian,
            4*max(scalar(0.25), max(mag(source.internalField())))
        );

        fvScalarMatrix PhiEqn
        (
            boostLaplacian*fvm::laplacian(Phi)
          - boostLaplacian*del2Phi
          + source
        );
        PhiEqn.setReference(0, scalar(0));
        solverPerformance sp = PhiEqn.solve();
        converged = sp.nIterations() == 0;
        
        // mimetic gradPhif for HessianVolumeRatio = 1
        if (HessianVolumeRatio >= 0.9)
        {
            gradPhif = TRiSK::divReconstruct(fvc::snGrad(Phi)*mesh.magSf());
            gradPhi = TRiSK::faceToCellMap(gradPhif);
        }
        // non-mimetic gradPhi if HessianVolumeRatio < 1
        else
        {
            gradPhi = fvc::reconstruct(fvc::snGrad(Phi)*mesh.magSf());
            gradPhif = fvc::interpolate(gradPhi);
            gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
                        *mesh.Sf()/mesh.magSf();
        }

        // create the new mesh depending of stencil size
        // Move the top and bottom patch point locations
        pointField newPoints = mesh.points();
        if (wideStencil)  // map grad phi from cell centres to edges
        {
            pointField gradPhiE = mte.interpolate(gradPhi);
            label halfNp = 0.5*newPoints.size();
            for(label ip = 0; ip < halfNp; ip++)
            {
                point np= expMapOnSphere
                (
                    0.5*(newPoints[ip] + newPoints[halfNp+ip]),
                    gradPhiE[ip]
                );
                np = unitVector(np);
                newPoints[ip] = np*mag(mesh.points()[ip]);
                newPoints[halfNp+ip] = np*mag(mesh.points()[halfNp+ip]);
            }
        }
        else // smaller stencil, recontructing gradPhi on points
        {
            // map gradPhi from cells onto points
            pointVectorField gradPhiP = fvc::faceToPointReconstruct(fvc::snGrad(Phi));

            label halfNp = 0.5*newPoints.size();
            for(label ip = 0; ip < halfNp; ip++)
            {
                point np= expMapOnSphere
                (
                    0.5*(newPoints[ip] + newPoints[halfNp+ip]),
                    0.5*(gradPhiP[ip] + gradPhiP[halfNp+ip])
                );
                np = unitVector(np);
                newPoints[ip] = np*mag(mesh.points()[ip]);
                newPoints[halfNp+ip] = np*mag(mesh.points()[halfNp+ip]);
            }
        }

        rMesh.movePoints(newPoints);
        
        // create the Voronoi generating points
        Vpoints = mesh.C() + gradPhi;
        
        // calculate the determinant of the Hessian
        Hessian = fvc::grad(gradPhif);
        forAll(detHess, cellI)
        {
            detHess[cellI] = det(diagTensor::one + Hessian[cellI]);
        }
//        volScalarField magGradPhi = mag(fvc::grad(Phi))/mag(mesh.C());
//        forAll(detHess, cellI)
//        {
//            detHess[cellI] = Foam::sin(magGradPhi[cellI])
//                             /mag(magGradPhi[cellI]);
//        }
        volRatio.internalField() =rMesh.V()/mesh.V();
        volRatio.correctBoundaryConditions();
        
        // combine the Hessian and the volume ratio
        detHess = (1-HessianVolumeRatio)*detHess + HessianVolumeRatio*volRatio;

        // Calculate the laplacian
        del2Phi = fvc::laplacian(Phi);

        // map to or calculate the monitor function to the new mesh and smooth
        monitorR = monitorFunc().map(rMesh, monitor);
        monitorNew.internalField() = monitorR.internalField();
        monitorNew.correctBoundaryConditions();

        // mean equidistribution
        equiDistMean = fvc::domainIntegrate(detHess)
                       /fvc::domainIntegrate(1/monitorNew);

        // Calculate the equidistribution for post processing
        equidist = detHess - equiDistMean/monitorNew;
        equidistVol = volRatio - equiDistMean/monitorNew;
        
        // Smooth the monitor function for the source term
        //monitorNew += 0.25*fvc::laplacian(sqr(1/mesh.deltaCoeffs()), monitorNew);
        
        // Calculate the source terms for the phi equation
        equiDistMean = fvc::domainIntegrate(detHess)
                       /fvc::domainIntegrate(1/monitorNew);
        source = detHess - equiDistMean/monitorNew;

        Info << "Time = " << runTime.timeName()
             << " determinant goes from " << min(detHess).value()
             << " to " << max(detHess).value()
             << " equidist goes from "
             << min(equidist.internalField())
             << " to " << max(equidist.internalField())
             << " cell volumes go from " << min(rMesh.V()).value() << " to "
             << max(rMesh.V()).value() << " ratio = "
             << (max(rMesh.V())/min(rMesh.V())).value()
             << " boostLaplacian = " << boostLaplacian << endl;

        if (converged)
        {
            runTime.writeAndEnd();
        }
        runTime.write();
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
