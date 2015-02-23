/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "snGradScheme.H"
#include "gaussLaplacianScheme.H"

namespace Foam
{
    namespace fvc
    {

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconstructd
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    // tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    // (
    //     new GeometricField<GradType, fvPatchField, volMesh>
    //     (
    //         IOobject
    //         (
    //             "volIntegrate("+ssf.name()+')',
    //             ssf.instance(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         inv(surfaceSum(sqr(mesh.Sf())/mesh.magSf()))
    //       & surfaceSum((mesh.Sf()/mesh.magSf())*ssf),
    //         zeroGradientFvPatchField<GradType>::typeName
    //     )
    // );

    surfaceVectorField dHat(mesh.delta()/mag(mesh.delta()));

    tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(surfaceSum(mesh.magSf()*sqr(dHat)))
          & surfaceSum(mesh.magSf()*dHat*ssf),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    treconField().correctBoundaryConditions();

    return treconField;
}
}
}


tmp<surfaceScalarField> nonOrthogonalCorrection
(
    const surfaceScalarField& snonCorrf
)
{
    const fvMesh& mesh = snonCorrf.mesh();

    surfaceScalarField nonCorrf
    (
        (mesh.deltaCoeffs()/mesh.nonOrthDeltaCoeffs())*snonCorrf
    );

    surfaceVectorField dHat(mesh.delta()/mag(mesh.delta()));

    surfaceVectorField nonCorrVec
    (
        fvc::interpolate(fvc::reconstructd<scalar>(nonCorrf))
    );

    nonCorrVec += dHat*(nonCorrf - (dHat&nonCorrVec));

    return (mesh.nonOrthCorrectionVectors() & nonCorrVec);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            twoPhaseProperties.correct();

            #include "alphaEqnSubCycle.H"
            interface.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
