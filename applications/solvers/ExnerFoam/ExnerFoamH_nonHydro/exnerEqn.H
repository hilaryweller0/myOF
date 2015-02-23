{
    Psi = pow(rho,(2*kappa-1)/(kappa-1))*pow(R/pRef*theta, kappa/(kappa-1));
    rho = Exner*Psi;
    rhof = fvc::interpolate(rho);

    surfaceScalarField G("G", 1/(1+offCentre*dt*muSponge));

    for (int corr=0; corr<nCorr; corr++)
    {
        // Explicit part of rate of change of flux
        dVdt = -H.delta() &
        (
            fvc::interpolate(fvc::div(U,u)) 
         + (2*rhof*(Omega^(Uf-Ug)))
        );

        // Vu is the flux without g or dpdz (to get the BCs right)
        Vu = G*(V.oldTime() + offCentre*dt*dVdt);

        // add the gravity terms and the hydrostatic pressure gradient
        surfaceScalarField gradPcoeff = G*offCentre*Cp*dt
                                        *rhof*thetaf*H.magd();
        V = Vu + G*offCentre*dt*rhof*gd - gradPcoeff*fvc::snGrad(ExnerH);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            // Flux in the normal direction (including off diagonal part of
            // pressure gradient)
            U = H.ddirToFlux(V)
              - H.ddirToFluxOffDiag(gradPcoeff*fvc::snGrad(ExnerNH));
 
            fvScalarMatrix ExnerEqn
            (
                fvm::ddt(Psi, ExnerNH)
              + fvc::ddt(Psi, ExnerH)
              + (1-offCentre)*divU
              + offCentre*fvc::div(U)
              - fvm::laplacian
                (
                    offCentre*H.Hdiag()*gradPcoeff/mesh.magSf(),
                    ExnerNH
                )
            );

            if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
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
                V -= gradPcoeff*fvc::snGrad(ExnerH);
                u = H.reconstructd(V/rhof);
                Uf = fvc::interpolate(u);
            }
        }
    }
    
    Exner == ExnerH + ExnerNH;
}

