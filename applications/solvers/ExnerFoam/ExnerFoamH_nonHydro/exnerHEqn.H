{
    // Solve hydrostatic balance for ExnerH
    fvScalarMatrix ExnerHEqn
    (
        fvc::div(H.Hdiag()*gd)
      - fvm::laplacian
        (
            H.Hdiag()*Cp*thetaf*H.magd()/mesh.magSf(),
            ExnerH
        )
    );
    ExnerHEqn.solve();
}
