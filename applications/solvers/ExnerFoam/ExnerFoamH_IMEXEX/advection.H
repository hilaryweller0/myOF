// SSP RK time-stepping as defined in GST01 (so far without sub-stepping)
{
    Info << "Advection RK steps" << flush;

    FixedList<scalar,3> alpha;
    alpha[0] = 1;
    alpha[1] = 0.25;
    alpha[2] = 2/3.;
    
    surfaceScalarField Umid = (1-offCentre)*U.oldTime() + offCentre*U;
    volScalarField divU = fvc::div(Umid);
    volScalarField rhoTheta = rho*theta;
    rhoTheta.oldTime();
    
    for(label irk = 0; irk < alpha.size(); irk++)
    {
        Info << ' ' << irk << flush;
        V = (1-alpha[irk])*V.oldTime() + alpha[irk]*V
          - alpha[irk]*dt*(fvc::interpolate(fvc::div(Umid,u)) & H.delta());
        rho = (1-alpha[irk])*rho.oldTime() + alpha[irk]*rho
            - alpha[irk]*dt*divU;
        rhoTheta = (1-alpha[irk])*rhoTheta.oldTime() + alpha[irk]*rhoTheta
                 - alpha[irk]*dt*fvc::div(Umid,theta);
        
        theta = rhoTheta/rho;
        u = H.reconstructd(V)/rho;
    }
    Info << endl;
}
