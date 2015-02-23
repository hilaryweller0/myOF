/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar using 3rd order Runge-Kutta

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // RK3 time-stepping
        T = T.oldTime() - dt/3*fvc::div(phi, T);
        T = T.oldTime() - dt/2*fvc::div(phi, T);
        T = T.oldTime() - dt*fvc::div(phi, T);

//        // RK2
//        T = T.oldTime() - dt*fvc::div(phi, T);
//        T = T.oldTime() - 0.5*dt*(fvc::div(phi, T.oldTime()) +  fvc::div(phi, T));
//        T = T.oldTime() - 0.5*dt*(fvc::div(phi, T.oldTime()) +  fvc::div(phi, T));

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
