// Sub step for advecting theta using and RK3 scheme
{
    dimensionedScalar dtsub = dt/nAdvectionSubSteps;

    Info << "Theta advection sub steps" << flush;
    
    rho = rho.oldTime();
    volScalarField rhoTheta(rho*theta);
    
    for(int subStep = 0; subStep < nAdvectionSubSteps; subStep++)
    {
        Info << " " << subStep << flush;
        
        divUt1 = fvc::div(U,theta);
        rhoTheta -= dtsub*divUt1;
        rho -= dtsub*divU;
        theta = rhoTheta/rho;
        
        divUt2 = fvc::div(U,theta);
        rhoTheta -= dtsub*(-0.75*divUt1 + 0.25*divUt2);
        rho += dtsub*0.5*divU;
        theta = rhoTheta/rho;
        
        divUt3 = fvc::div(U,theta);
        rhoTheta -= dtsub*(-1/12.*(divUt1 + divUt2) + 2/3.*divUt3);
        rho -= dtsub*0.5*divU;
        theta = rhoTheta/rho;
    }
    
    Info << endl;
}
