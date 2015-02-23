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
    setBalancedJet

Description
    Sets initial conditions for shallowWaterFoamRK, SWE solver using TRiSK.
    Initial conditions for steady state, barotropically unstable jet

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    const Switch overRideCellGeometry
    (
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "overRideCellGeometry", Switch(false)
        )
    );
    const Switch sampleh
    (
        runTime.controlDict().lookupOrDefault<Switch>("sampleh", Switch(false))
    );

    Foam::fvMeshWithDual mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        overRideCellGeometry,
        fvMeshWithDual::SPHERICALDIST
    );
    
    Foam::Info
        << "Create dualMesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMeshWithDual dualMesh
    (
        Foam::IOobject
        (
            "dualMesh",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        mesh,
        overRideCellGeometry
    );
    mesh.setDual(dualMesh);

    TRiSKData trisk(mesh);
    #define HDIAG trisk.circToFluxDiag()

    const Switch conserveEnergy
    (
        mesh.solutionDict().lookup("conserveEnergy")
    );

    Info << "\nReading earthProperties" << endl;

    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "Reading properties for initial conditions" << endl;
    const dimensionedScalar& magg = mesh.magg();
    const dimensionedVector Omega = mesh.Omega();
    const dimensionedScalar Rsphere("Rsphere", dimLength, mag(mesh.C()[0]));
    const scalar degToRad = constant::mathematical::pi/180.;
    const dimensionedScalar h0(earthProperties.lookup("h0"));
    const dimensionedScalar umax(earthProperties.lookup("umax"));
    const scalar lat0 = readScalar(earthProperties.lookup("lat0"))*degToRad;
    const scalar lat1 = readScalar(earthProperties.lookup("lat1"))*degToRad;
    const scalar dlat = readScalar(earthProperties.lookup("dlat"))*degToRad;

    const scalar en = Foam::exp(-4/sqr(lat1 - lat0));
    const scalar aue = Rsphere.value()*umax.value()/en;

    const volScalarField f
    (
        IOobject("f", runTime.constant(), dualMesh),
        2*(Omega & dualMesh.C())/mag(dualMesh.C())
    );

    Info << "Creating ff, Coriolis on the primal face centres" << endl;
    surfaceScalarField ff
    (
        IOobject("ff", runTime.timeName(), mesh),
        2*(Omega & mesh.Cf())/mag(mesh.Cf())
    );
    ff.write();
    
    Info << "Creating ff, Coriolis on the dual face centres" << endl;
    surfaceScalarField ffd
    (
        IOobject("ff", runTime.timeName(), dualMesh),
        mesh.dualMap(ff)
    );
    ffd.write();

    Info << "Creating field streamFunc by numerical integration\n" << endl;
    volScalarField streamFunc
    (
        IOobject("streamFunc", runTime.timeName(), dualMesh),
        dualMesh,
        dimensionedScalar("streamFunc", dimensionSet(0,2,-1,0,0), scalar(0))
    );
    forAll(streamFunc, ip)
    {
        const scalar lat = Foam::asin(dualMesh.rHat()[ip].z());

        label nlat = ceil(mag(lat)/dlat);
        scalar dlatp = lat/max(scalar(1), scalar(nlat));
        for(label il = 0; il < nlat; il++)
        {
            const scalar latp = dlatp*(il+0.5);
            if (latp > lat0 - SMALL && latp < lat1 + SMALL)
            {
                streamFunc[ip] += Foam::exp(1/((latp - lat0)*(latp - lat1)));
            }
        }
        streamFunc[ip] *= -dlatp*aue;
    }

    Info << "Creating vS by integrating stream function\n" << endl;
    surfaceScalarField vS
    (
        IOobject("vS", runTime.timeName(), dualMesh),
        fvc::snGrad(streamFunc)*dualMesh.magSf()
    );

    // Velocity vector field
    surfaceVectorField Ufd
    (
        "Ufd",
        TRiSK::faceReconstruct(vS/dualMesh.magSf(), dualMesh.jdir())
    );
    Ufd -= (Ufd & dualMesh.rHatf())*dualMesh.rHatf();
    surfaceVectorField Uf("Uf", dualMesh.dualMap(Ufd));

    Info << "Creating h to be in discrete geostrophic balance with vS\n" << endl;
    volScalarField h
    (
        IOobject("h", runTime.timeName(), mesh),
        mesh,
        h0
    );
    
    // initial analytic h profile if sampleh
    if (sampleh)
    {
        forAll(h, cellI)
        {
            h[cellI] = 0;
            scalar lat = mesh.lat()[cellI];
            label nlat = ceil(mag(lat)/dlat);
            scalar dlatp = lat/max(scalar(1), scalar(nlat));
            for (label i = 0; i < nlat; i++)
            {
                scalar latp = dlatp*(i+0.5);
                scalar u = 0;
                if (latp > lat0-SMALL && latp < lat1+SMALL)
                {
                    u = (umax/en * Foam::exp(1/((latp - lat0)*(latp - lat1)))).value();
                }
                scalar fp = 2*Omega[2].value()*Foam::sin(latp);
                h[cellI] -= u*(Rsphere.value()*fp + Foam::tan(latp)*u);
            }
            h[cellI] *= dlatp/magg.value();
        }
        // set the level of the height so that the global mean is h0
        dimensionedScalar meanh = sum(mesh.V()*h)/sum(mesh.V());
        h += h0 - meanh;
    }
    else // find discretely balanced h
    {
        bool converged = false;

        volScalarField KE = TRiSK::KE(TRiSK::circToFlux(vS), vS);
        volScalarField pv("pv", (TRiSK::curl(vS)+ f)/TRiSK::primalToDualCellMap(h));
        surfaceScalarField hf = fvc::interpolate(h);
        surfaceScalarField hfd = mesh.dualMap(hf);
        surfaceScalarField phi = TRiSK::circToFlux(hfd*vS);
        surfaceScalarField phid = Ufd & dualMesh.Sf();
        surfaceScalarField pvf = fvc::interpolate(pv);
        surfaceScalarField pvFlux = pvf*TRiSK::perp(phi);
        if (conserveEnergy) pvFlux = 0.5*
        (
            pvFlux
          + TRiSK::perp(dualMesh.dualMap(pvf)*phi)
        );

        for(label it = 0; it < 40 && !converged; it++)
        {
            Info << "\nIteration " << it << endl;
        
            // set the level of the height so that the global mean is h0
            dimensionedScalar meanh = fvc::domainIntegrate(h)/sum(mesh.V());
            Info << "meanh - h0 = " << meanh.value() - h0.value() << endl;
            h += h0 - meanh;

            // Solve h equation to make h in geostrophic balance with un
            fvScalarMatrix hEqn
            (
                fvm::laplacian(magg*HDIAG, h)
             == -fvc::div
                 (
                    TRiSK::circToFlux
                    (
                        pvFlux - mesh.dualFluxMap(fvc::snGrad(KE)*mesh.magSf())
                    )
                  - magg*TRiSK::circToFluxOffDiag
                    (
                        mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf())
                    )
                 )
            );
            //hEqn.setReference(0, h[0]);
            
            converged = hEqn.solve(mesh.solver(h.name() + "Final")).nIterations()
                          <= 1;

            pv = (TRiSK::curl(vS)+ f)/TRiSK::primalToDualCellMap(h);
            hf = fvc::interpolate(h);
            hfd = mesh.dualMap(hf);
            phi = TRiSK::circToFlux(hfd*vS);
            phid = Ufd & dualMesh.Sf();
            pvf = fvc::interpolate(pv);
            pvFlux = pvf*TRiSK::perp(phi);
            if (conserveEnergy)
            {
                pvFlux = 0.5*
                (
                    pvFlux
                  + TRiSK::perp(dualMesh.dualMap(pvf)*phi)
                );
            }
        }
    }
    h.write();
    Uf.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
