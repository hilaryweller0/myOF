// Sub step advection using and RK3 scheme
{
    dimensionedScalar dtsub = dt/nAdvectionSubSteps;

    // mid-time-step value for rho (makes it unstable?)
    rho = rho.oldTime() - (1-offCentre)*dt*divPhi;
    rhof = fvc::interpolate(rho);

    Info << "Advection sub steps" << flush;
    for(int subStep = 0; subStep < nAdvectionSubSteps; subStep++)
    {
        Info << " " << subStep << flush;
        U = fvc::reconstruct(phi.oldTime()/rhof);
        divUu1 = fvc::interpolate(fvc::div(phi,U)) & mesh.Sf();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu1 -= diffusionCoeff*rhof*
                      (mesh.Sf() & fvc::interpolate(fvc::laplacian(U)));
        }
        phi.oldTime() -= dtsub*divUu1;
        
        U = fvc::reconstruct(phi.oldTime()/rhof);
        divUu2 = fvc::interpolate(fvc::div(phi,U)) & mesh.Sf();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu2 -= diffusionCoeff*rhof*
                      (mesh.Sf() & fvc::interpolate(fvc::laplacian(U)));
        }
        phi.oldTime() += dtsub*(0.75*divUu1 - 0.25*divUu2);
        
        U = fvc::reconstruct(phi.oldTime()/rhof);
        divUu3 = fvc::interpolate(fvc::div(phi,U)) & mesh.Sf();
        if (mag(diffusionCoeff.value()) > SMALL)
        {
            divUu3 -= diffusionCoeff*rhof*
                      (mesh.Sf() & fvc::interpolate(fvc::laplacian(U)));
        }
        phi.oldTime() += dtsub*(1/12.*(divUu1 + divUu2) - 2/3.*divUu3);
    }
    
    // update velocity variables
    U = fvc::reconstruct(phi.oldTime()/rhof);
    Uf = fvc::interpolate(U);
    Info << endl;
}
