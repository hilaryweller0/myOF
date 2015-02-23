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

    surfaceScalarField sndKg
    (
        fvc::interpolate(interface.sigmaK())
       *fv::snGradScheme<scalar>::sndGrad(alpha1)
      - ghf*fv::snGradScheme<scalar>::sndGrad(rho)
    );

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phiCorr
      + rAUf
       *(
            fvc::interpolate(rho.oldTime())*phi.oldTime()/runTime.deltaT()
          - (fvc::interpolate(UEqn & U) & mesh.Sf())
          + sndKg*mesh.magSf()
        )
    );

    phi = phiHbyA;

    while (pimple.correctNonOrthogonal())
    {
        surfaceScalarField nonOrthCorr
        (
            nonOrthogonalCorrection
            (
                rAUf*
                (
                    sndKg
                  - fv::snGradScheme<scalar>::sndGrad(p_rgh)
                )
            )*mesh.magSf()
        );

        fvScalarMatrix p_rghEqn
        (
            fv::gaussLaplacianScheme<scalar, scalar>::fvmLaplacianUncorrected
            (
                rAUf*mesh.magSf(),
                mesh.nonOrthDeltaCoeffs(),
                p_rgh
            )
            //fvm::laplacian(rAUf, p_rgh)
         == fvc::div(phiHbyA + nonOrthCorr)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + nonOrthCorr - p_rghEqn.flux();

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
