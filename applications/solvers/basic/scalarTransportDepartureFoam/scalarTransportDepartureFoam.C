/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar using a forwrd in time discretisation
    based on a single departure point per face

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToPointField.H"
#include "rbfFit.H"
#include "polyFit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"
    
    // Departure points for each face
    surfaceVectorField departurePoints
    (
        IOobject("departurePoints", runTime.timeName(), mesh),
        mesh.Cf() - 0.5*Uf*dt
    );
    departurePoints.write();
    
    // Class for interpolating onto departure points
    meshToPointField<polyFit<oONE> > meshToDep
    (
        departurePoints, mesh, meshToPointField<polyFit<oONE> >::CELLCELLS
    );
    
    // Interpolate the velocity onto the departure points, calculate the flux and
    // make it divergence free
    //Uf.internalField() = meshToDep.interpolate(U);
    phi = Uf & mesh.Sf();

    for(int icorr = 0; icorr < 12; icorr++)
    {
        fvScalarMatrix pEqn(fvm::laplacian(p) + fvc::div(phi));
        pEqn.setReference(0, scalar(0));
        pEqn.solve();
        if (icorr == 11) phi += pEqn.flux();
    }
    
    //Uf = linearInterpolate(fvc::reconstruct(phi));
    Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
    Uf.write();
    
//    volVectorField TU("TU", T*U);
//    surfaceVectorField TUf("TUf", linearInterpolate(TU));

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << flush;

        Tf.internalField() = meshToDep.interpolate(T);
//        TU = T*U;
//        TUf.internalField() = meshToDep.interpolate(TU);

        T = T.oldTime() - dt*fvc::div(phi*Tf);
        
        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

