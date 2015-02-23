{
    surfaceScalarField phiCorrCoeff
    (
        0.1 // MUST be <= 1
       *mag((fvc::interpolate(U) & mesh.Sf()) - phi)
       /(
            mag(fvc::interpolate(U) & mesh.Sf())
          + mag(phi)
          + dimensionedScalar("0", phi.dimensions(), 1e-15)
        )
    );

    Info<< "phiCorrCoeff " << phiCorrCoeff.average().value() << endl;

    surfaceScalarField phiCorr
    (
        phiCorrCoeff*((fvc::interpolate(U) & mesh.Sf()) - phi)
    );

    surfaceScalarField rAUf
    (
        "Dp",
        runTime.deltaT()/fvc::interpolate(rho)
    );

    surfaceScalarField phigK
    (
        rAUf
       *(
           fvc::interpolate(interface.sigmaK())*fv::snGradScheme<scalar>::sndGrad(alpha1)
         - ghf*fv::snGradScheme<scalar>::sndGrad(rho)
        )
    );

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phiCorr
      + rAUf
       *(
            fvc::interpolate(rho.oldTime())*phi.oldTime()/runTime.deltaT()
          - (fvc::interpolate(UEqn & U) & mesh.Sf())
        )
      + phigK*mesh.magSf()
    );

    phi = phiHbyA;

    while (pimple.correctNonOrthogonal())
    {
        surfaceScalarField nonCorrf
        (
            phigK
          - rAUf*
            (
                fv::snGradScheme<scalar>::sndGrad(p_rgh)
            )
        );

        surfaceVectorField nonCorrVec
        (
            fvc::interpolate(fvc::reconstructd<scalar>(nonCorrf))
        );

        nonCorrVec += dHat*(nonCorrf - (dHat&nonCorrVec));

        surfaceScalarField nonCorr
        (
            (mesh.nonOrthCorrectionVectors() & nonCorrVec)*mesh.magSf()
        );

        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA + nonCorr)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + nonCorr - p_rghEqn.flux();

            U = fvc::reconstruct(phi);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
