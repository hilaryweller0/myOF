/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    vorticityDual

Description
    Calculates the vorticity on the dual mesh

\*---------------------------------------------------------------------------*/

#include "fvMeshWithDual.H"
#include "TRiSK.H"
#include "fvCFD.H"

using namespace Foam;
//using namespace mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createMeshWithDual.H"
    #include "createDualMesh.H"
    word Htype
    (
        mesh.solutionDict().lookupOrDefault<word>("nonOrthogHtype", "diagonal")
    );

    const volScalarField f
    (
        IOobject("f", runTime.constant(), dualMesh),
        2*(dualMesh.Omega() & dualMesh.C())/mag(dualMesh.C())
    );
    
    Info << "mesh.Omega = " << mesh.Omega() << " dualMesh.Omega = " << dualMesh.Omega() << endl;
    
    const volScalarField fprimal
    (
        IOobject("f", runTime.constant(), mesh),
        2*(mesh.Omega() & mesh.C())/mag(mesh.C())
    );
    
    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "Reading field pv\n" << endl;
        volScalarField pv
        (
            IOobject("pv", runTime.timeName(), dualMesh, IOobject::MUST_READ),
            dualMesh
        );
        
        Info << "Reading field h\n" << endl;
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
            mesh,
            dimensionedScalar("H", dimless, scalar(1))
        );
        
        Info<< "Reading field Uf\n" << endl;
        surfaceVectorField Uf
        (
            IOobject("Uf", runTime.timeName(), mesh, IOobject::MUST_READ),
            mesh
        );
        
        surfaceVectorField Ufd("Ufd", mesh.dualMap(Uf));
        
        // The velocity flux in the dual perpendicular direction (the d direction)
        surfaceScalarField vS
        (
            IOobject
            (
                "vS", runTime.timeName(), dualMesh,
                IOobject::READ_IF_PRESENT
            ),
            (Ufd & dualMesh.jdir())*dualMesh.magSf()
        );
        if (Htype == "diagonal")
        {
            vS = -(Ufd & mesh.signedDualMap(mesh.idir()))*dualMesh.magSf();
        }

        Info << "Vorticity on the dual mesh\n" << endl;
        volScalarField hv("h",TRiSK::primalToDualCellMap(h));
        volScalarField vorticity("vorticity", pv*hv - f);
        vorticity.write();
        
        Info << "Vorticity on the primal mesh\n" << endl;
        volScalarField vorticityp("vorticity", TRiSK::primalToDualCellMap(vorticity));
        vorticityp.write();
        
        Info << "Divergence on primal and dual meshes\n" << endl;
//        volScalarField divu("divu", fvc::div(un*mesh.magSf()));
//        divu.write();
        
        surfaceScalarField phi
        (
            TRiSK::circToFlux(mesh.dualMap(fvc::interpolate(h))*vS)
        );
            
        volScalarField divhu("divhu", fvc::div(phi));
        divhu.write();

        volScalarField divhud("divhu", fvc::div(TRiSK::perp(phi)));
        divhud.write();

//        // KE on the primal mesh
//        volScalarField KE("KE", TRiSK::ke(un));
//        KE.write();
        
//        // enstrophy on the dual mesh
//        hv.write();
//        volScalarField enstrophy("enstrophy", sqr(vorticity+f)/hv);
//        enstrophy.write();
    }
}


// ************************************************************************* //

