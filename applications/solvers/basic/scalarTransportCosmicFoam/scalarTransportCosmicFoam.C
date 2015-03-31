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
    Solves a transport equation for a passive scalar using COSMIC dimension
    splitting

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
    #include "create1dStencils.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int corr = 0; corr < 2; corr++)
        {
            // 1d Advective updates for cross terms
            gradT = fvc::grad(fvc::interpolate(T));
            volScalarField Tax = T.oldTime() - dt*u*gradT.component(vector::X);
            //volScalarField Tay = T.oldTime() - dt*v*gradT.component(vector::Y);
            volScalarField Taz = T.oldTime() - dt*w*gradT.component(vector::Z);
        
            // Conservative updates using cross terms
            volVectorField graduT = fvc::grad(uf*fvc::interpolate(T+Taz, "interpolate(T)"));
            volVectorField gradwT = fvc::grad(wf*fvc::interpolate(T+Tax, "interpolate(T)"));
            T = T.oldTime() - 0.5*dt*(graduT.component(vector::X))
                            - 0.5*dt*(gradwT.component(vector::Z));
        }

        T.correctBoundaryConditions();
        
        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
