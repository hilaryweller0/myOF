{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    thetaf = fvc::interpolate(theta);

    // Other Explicit part of rate of change of flux
    dVdt = -H.delta() & (2*rhof*(Omega^(Uf-Ug)));

    // Vu is the flux without g or dpdz (to get the BCs right)
    Vu = G*(V.oldTime() + offCentre*dt*dVdt);

    // add the gravity terms
    V = Vu + G*offCentre*dt*rhof*gd;

    // Flux in the normal direction (including off diagonal part of pressure
    // gradient)
    gradPcoeff = G*offCentre*Cp*dt*rhof*thetaf*H.magd();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        U = H.ddirToFlux(V)
          - H.ddirToFluxOffDiag(gradPcoeff*fvc::snGrad(Exner));

        fvScalarMatrix ExnerEqn
        (
            fvm::ddt(Psi, Exner)
          + (1-offCentre)*divU
          + offCentre*fvc::div(U)
          - fvm::laplacian(offCentre*H.Hdiag()*gradPcoeff/mesh.magSf(), Exner)
        );

        if(nonOrth == nNonOrthCorr && icorr == nCorr-1)
        {
            ExnerEqn.solve(mesh.solver(Exner.name() + "Final"));
        }
        else
        {
            ExnerEqn.solve(mesh.solver(Exner.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            U += ExnerEqn.flux()/offCentre;
            //V -= gradPcoeff*fvc::snGrad(Exner);
            u = fvc::reconstruct(U)/rho;
            Uf = fvc::interpolate(u);
        }
    }
}
