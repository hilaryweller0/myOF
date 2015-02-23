/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    setAtmosProfile

Description
    Set Exner and theta for constant Brunt Vaisalla frequency

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ExnerTheta.H"
#include "ITstream.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readEnvironmentalProperties.H"
#   include "readThermoProperties.H"
#   include "readInitialProperties.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Setting theta_init\n" << endl;
    if (mag(g.value()) > SMALL)
    {
        theta_init == T0*exp(-sqr(BruntV)*(mesh.C() & ghat)/mag(g));
    }
    theta_init.write();

    Info<< "Setting hydrostatically balanced Exner_init\n" << endl;
    if (mag(BruntV.value()) > SMALL)
    {
        Exner_init == 1 - kappa*magSqr(g)
                      /(T0*R*sqr(BruntV))*(1 - T0/theta_init);
    }
    else
    {
        Exner_init == 1 + kappa*(mesh.C() & g)/(theta_init*R);
    }
    Exner_init.write();

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
