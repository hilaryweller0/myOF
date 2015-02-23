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
    sphere. Reads in a monitor function to define the mesh density.

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
#include "VectorSpaceFunctions.H"
//#include "meshToMesh.H"
//#include "volPointInterpolation.H"
#include "monitorFunction.H"
#include "meshToPoints.H"
#include "meshToEdges.H"
#include "polyFit.H"
//#include "VoronoiSphereMeshing.H"
//#include "InitialPointsRaw.H"

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
//    fvMesh mesh
//    (
//        Foam::IOobject
//        (
//            meshRegion, runTime.timeName(), runTime, IOobject::MUST_READ
//        )
//    );
    
    // Create the mesh to be moved
    fvMesh rMesh
    (
        Foam::IOobject
        (
            "rMesh", runTime.timeName(), runTime,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        )
    );

    // maps for mapping from cell centres to points
    //meshToPoints<polyFit<1> > mtp(mesh);
    meshToEdges<polyFit<1> > mte(mesh);
    
    // Open control dictionary
    IOdictionary controlDict
    (
        IOobject
        (
            args.executable() + "Dict", runTime.system(), runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
//    // Dictionary for the Voronoi meshing
//    IOdictionary Vdict
//    (
//        IOobject
//        (
//            "VoronoiSphereMeshDict", runTime.system(), runTime,
//            IOobject::MUST_READ_IF_MODIFIED
//        )
//    );
    
    // The monitor funciton
    autoPtr<monitorFunction> monitorFunc(monitorFunction::New(controlDict));
    
    // Read the ratio of the Hessian to volume to use
    const scalar HessianVolumeRatio
    (
        controlDict.lookupOrDefault<scalar>("HessianVolumeRatio", scalar(0))
    );
    
    #include "createFields.H"
    //#include "nFacesPerPoint.H"

//    // create the initial Voronoi surface mesh
//    InitialPointsRaw ivPoints(mesh.C());
//    VoronoiSphereMeshing vMeshing(Vdict, ivPoints);
//    const faceList dualFace = vMeshing.dualFaces();
//    PrimitivePatch<face, List, pointField> Vpatch0
//    (
//        dualFace, vMeshing.dualPoints(dualFace)
//    );
//    scalarList VpatchFaceAreas0(Vpatch0.size());
//    forAll(VpatchFaceAreas0, faceI)
//    {
//        VpatchFaceAreas0[faceI] = Vpatch0[faceI].mag(Vpatch0.points());
//    }

    // Diffusivity to control stability
    //dimensionedScalar diffusivity("diffusivity", dimensionSet(0,2,-1,0,0),scalar(0));
    dimensionedScalar diffusivity(controlDict.lookup("diffusivity"));
    dimensionedScalar rdiffusivity(controlDict.lookup("rdiffusivity"));

    // Store 4 initial residuals in a row to check for non-convergence
    scalarList initResids(4, GREAT);

    // Use time-steps instead of iterations to solve the Monge-Ampere eqn
    bool converged = false;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << flush << nl;
        
        equiDistMean = fvc::domainIntegrate(detHess)
                      /fvc::domainIntegrate(1./monitorNew);

        fvScalarMatrix PhiEqn
        (
            rdiffusivity*fvm::ddt(Phi)
          - (1+diffusivity)*fvm::laplacian(Phi)
          + equiDistMean/monitorNew
          + del2Phi
          - detHess
        );
        //PhiEqn.setReference(0, scalar(0));
        solverPerformance sp = PhiEqn.solve();
        converged = sp.nIterations() == 0;
        
        // Check that initial residuals are not diverging
        for(label ir = 0; ir < initResids.size()-1; ir++)
        {
            initResids[ir] = initResids[ir+1];
        }
        initResids.last() = sp.initialResidual();
        int diverging = 1;
        for(label ir = 0; ir < initResids.size()-1; ir++)
        {
            diverging = diverging * int(initResids.last() > initResids[ir]);
        }
        if (diverging) converged = true;

        // mimetic gradPhif for HessianVolumeRatio = 1
        if (HessianVolumeRatio >= 1-SMALL)
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

        // create the new mesh
        // map gradPhi from cells onto internal edges
        pointField gradPhiE = mte.interpolate(gradPhi);
        //pointField newPoints = mesh.points() + mtp.interpolate(gradPhi);
        
//        // map gradPhi from faces onto points
//        pointField newPoints = mesh.points();
//        forAll(gradPhif, faceI)
//        {
//            const face& f = mesh.faces()[faceI];
//            forAll(f, fi)
//            {
//                newPoints[f[fi]] += gradPhif[faceI]/scalar(nFacesPerPoint[f[fi]]);
//            }
//        }
        
        // Move the top and bottom patch point locations
        pointField newPoints = mesh.points();
        label halfNp = 0.5*newPoints.size();
        for(label ip = 0; ip < halfNp; ip++)
        {
            point np = 0.5*(newPoints[ip] + newPoints[halfNp+ip]) + gradPhiE[ip];
            np = unitVector(np);
            newPoints[ip] = np*mag(mesh.points()[ip]);
            newPoints[halfNp+ip] = np*mag(mesh.points()[halfNp+ip]);
//            newPoints[ip] = 0.5*(newPoints[ip] + newPoints[halfNp+ip]);
//            newPoints[halfNp+ip] = newPoints[ip];
        }
//        newPoints *= mag(mesh.points())/mag(newPoints);
        
        rMesh.movePoints(newPoints);
        
        // create the Voronoi generating points
        Vpoints = mesh.C() + gradPhi;
        
//        // create the Voronoi surface mesh
//        InitialPointsRaw ivPoints(Vpoints);
//        VoronoiSphereMeshing vMeshing(Vdict, ivPoints);
//        const faceList dualFace = vMeshing.dualFaces();
//        PrimitivePatch<face, List, pointField> Vpatch
//        (
//            dualFace, vMeshing.dualPoints(dualFace)
//        );
//        scalarList VpatchFaceAreas(Vpatch.size());
//        forAll(VpatchFaceAreas, faceI)
//        {
//            VpatchFaceAreas[faceI] = Vpatch[faceI].mag(Vpatch.points());
//        }

        // calculate the determinant of the Hessian
        Hessian = fvc::grad(gradPhif);
        Hessian = 0.5*(Hessian + Hessian.T());
        del2Phi = fvc::laplacian(Phi);
        forAll(detHess, cellI)
        {
            detHess[cellI] = det(diagTensor::one + Hessian[cellI]);
        }
        volRatio.internalField() =rMesh.V()/mesh.V();
        volRatio.correctBoundaryConditions();
        
        detHess = (1-HessianVolumeRatio)*detHess + HessianVolumeRatio*volRatio;
//        detHess.internalField() = VpatchFaceAreas/VpatchFaceAreas0;

        boostLaplacian = detHess/(1+del2Phi);

        // map to or calculate the monitor function to the new mesh
        //monitorR = monitorFunc().map(rMesh, monitor);
        monitorR.internalField() = monitorFunc().map(Vpoints, monitor);

        monitorNew.internalField() = monitorR.internalField();
        
        equidist = monitorNew*detHess;
        //equidist = monitorNew*volRatio;
        
        Info << "Time = " << runTime.timeName()
             << " determinant goes from " << min(detHess).value()
             << " to " << max(detHess).value()
             << " equidist goes from "
             << min(equidist.internalField())
             << " to " << max(equidist.internalField())
             << " cell volume ratios go from "
             << min(volRatio).value() << " to "
             << max(volRatio).value() << " ratio = "
             << max(volRatio).value()/min(volRatio).value() << endl;

        if (/*min(detHess).value() < 0 ||*/ converged)
        {
            runTime.writeAndEnd();
        }
        runTime.write();
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
