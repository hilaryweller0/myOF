// SSP RK time-stepping as defined in GST01 (with sub-stepping)
{
    Info << "Advection RK steps and sub steps" << endl;

    dimensionedScalar dtsub = dt/nAdvectionSubSteps;

    FixedList<scalar,3> alpha;
    alpha[0] = 1;
    alpha[1] = 0.25;
    alpha[2] = 2/3.;
    
    divU = fvc::div(U);
    volScalarField rhoTheta = rho*theta;
    rhoTheta.oldTime();
    Vu = V.oldTime();
    u = H.reconstructd(Vu)/rho;
    
    for(int subStep = 0; subStep < nAdvectionSubSteps; subStep++)
    {
        Info << subStep << ": RK step" << flush;

        for(label irk = 0; irk < alpha.size(); irk++)
        {
            Info << ' ' << irk << flush;
            Vu = (1-alpha[irk])*V.oldTime() + alpha[irk]*Vu
               - alpha[irk]*dtsub*(fvc::interpolate(fvc::div(U,u))&H.delta());
            rho = (1-alpha[irk])*rho.oldTime() + alpha[irk]*rho
                - alpha[irk]*dtsub*divU;
            rhoTheta = (1-alpha[irk])*rhoTheta.oldTime() + alpha[irk]*rhoTheta
                     - alpha[irk]*dtsub*fvc::div(U,theta);
        
            theta = rhoTheta/rho;
            u = H.reconstructd(Vu)/rho;
        }
        V.oldTime() = Vu;
        rho.oldTime() = rho;
        rhoTheta.oldTime() = rhoTheta;
        Info << endl;
    }
}
