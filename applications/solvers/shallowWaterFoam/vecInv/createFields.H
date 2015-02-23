Info<< "Reading field Uf\n" << endl;
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
surfaceVectorField Ufd
(
    IOobject
    (
        "Ufd",
        runTime.timeName(),
        dualMesh
    ),
    mesh.dualMap(Uf)
);

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
if (trisk.Htype() == TRiSKData::DIAGONAL)
{
    vS = -(Ufd & mesh.signedDualMap(mesh.idir()))*dualMesh.magSf();
}
vS.oldTime();

Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field h0 if present\n" << endl;
const volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("h0", dimLength, scalar(0))
);
Info << "max h0 = " << max(h0) << endl;

const surfaceScalarField gGradh0 = magg
                             *mesh.dualFluxMap(fvc::snGrad(h0)*mesh.magSf());

// h on the faces and on the dual
surfaceScalarField hf
(
    IOobject("hf", runTime.timeName(), mesh),
    TRiSK::conserveInterp(h)
);
if (linear) hf == H;
volScalarField hv(IOobject("hv", runTime.timeName(), dualMesh), dualMesh, H);
if (!linear) hv = TRiSK::primalToDualCellMap(h);
surfaceScalarField hfd("hfd", mesh.dualMap(hf));

// The flux. Ie the velocity in the primal normal direction (hf*Sf)
surfaceScalarField phi
(
    IOobject
    (
        "phi", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    TRiSK::circToFlux(hfd*vS)
);
phi.oldTime();

// Flux on the dual mesh (for up-winded interpolation schemes)
surfaceScalarField phid
(
    IOobject("phid", runTime.timeName(), dualMesh),
    Ufd & dualMesh.Sf()
);

// recalculate hf, hfd and phi using the correct interpolation scheme
if (!linear)
{
    hf = fvc::interpolate(h);
    hfd = mesh.dualMap(hf);
    phi = TRiSK::circToFlux(hfd*vS);
}

const volScalarField f
(
    IOobject("f", runTime.constant(), dualMesh),
    2*(Omega & dualMesh.C())/mag(dualMesh.C())
);

const surfaceScalarField ff
(
    IOobject
    (
        "ff",
        runTime.timeName(),
        dualMesh,
        IOobject::READ_IF_PRESENT
    ),
    2*(Omega & dualMesh.Cf())/mag(dualMesh.Cf())
);

volScalarField vorticityByH
(
    IOobject("vorticityByH", runTime.timeName(), dualMesh),
    TRiSK::curl(vS)/hv
);

volScalarField pv
(
    IOobject
    (
        "pv", runTime.timeName(), dualMesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    vorticityByH + f/hv
//    (f + TRiSK::curl(vS))/hv
);
if (linear) pv = f/H;
pv.write();
surfaceScalarField pvf("pvf", ff/hfd);

if (!linear)
{
    if (interpolateVorticity)
    {
        pvf = ff/hfd;
        if (!linear)
        {
            pvf += fvc::interpolate(vorticityByH);
        }
    }
    else    
    {
        pvf = fvc::interpolate(pv);
    }
}

surfaceScalarField pvFlux
(
    IOobject("pvFlux", runTime.timeName(), dualMesh),
    pvf*TRiSK::perp(phi)
);
if (useTriskPerp && conserveEnergy)
{
    pvFlux = 0.5*(pvFlux + TRiSK::perp(dualMesh.dualMap(pvf)*phi));
}
else if (!useTriskPerp)
{
    pvFlux = pvf*TRiSK::fluxMap(phi, mesh.idir(), dualMesh.idir());
}

volScalarField KE
(
    IOobject
    (
        "KE", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    TRiSK::KE(TRiSK::circToFlux(vS), vS)
);
KE.write();

surfaceScalarField dvSdt
(
    IOobject("dvSdt", runTime.timeName(), dualMesh),
    -pvFlux
    + gGradh0
    + magg*mesh.dualFluxMap(fvc::snGrad(h)*mesh.magSf())
);
if (!linear) dvSdt += mesh.dualFluxMap(fvc::snGrad(KE)*mesh.magSf());
dvSdt.oldTime();

// Create and write out Uf consistent with the discretisation
Ufd = TRiSK::faceReconstruct
(
    vS/dualMesh.magSf(), dualMesh.jdir()
);
Ufd -= (Ufd & dualMesh.rHatf())*dualMesh.rHatf();
//Uf = dualMesh.dualMap(Ufd);
//{
//    surfaceVectorField Uf0("Uf0", Uf);
//    Uf0.write();
//}

surfaceScalarField del2u
(
    IOobject("del2u", runTime.timeName(), dualMesh),
  - mesh.signedDualMap
    (
        fvc::snGrad(fvc::div(TRiSK::circToFlux(vS)))
    )
  - fvc::snGrad(TRiSK::curl(vS))
);
surfaceVectorField del2U
(
    IOobject
    (
        "del2U", runTime.timeName(), dualMesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    TRiSK::faceReconstruct(del2u, dualMesh.jdir())
);
del2U.write();

