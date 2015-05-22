//V.oldTime() -= dt*
//(
//    H.delta() &
//    (
//        (1-offCentre)*fvc::interpolate(fvc::div(U.oldTime(),u.oldTime()))
//      + offCentre*fvc::interpolate(fvc::div(U,u))
//    )
//);


// SSP RK time-stepping as defined in GST01 (so far without sub-stepping)
{
    Info << "Advection RK steps" << flush;

    FixedList<scalar,3> alpha;
    alpha[0] = 1;
    alpha[1] = 0.25;
    alpha[2] = 2/3.;
    
    Vu = V.oldTime();
    u = H.reconstructd(Vu)/rho;
    
    for(label irk = 0; irk < alpha.size(); irk++)
    {
        Info << ' ' << irk << flush;
        Vu = (1-alpha[irk])*V.oldTime() + alpha[irk]*Vu
          - alpha[irk]*dt*(fvc::interpolate(fvc::div(U,u)) & H.delta());
        u = H.reconstructd(Vu)/rho;
    }
    V.oldTime() = Vu;
    Info << endl;
}
