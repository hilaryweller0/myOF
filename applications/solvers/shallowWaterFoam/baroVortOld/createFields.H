// The total velocity vector on the faces of the primal mesh
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

// The velocity in the primal normal direction (*Sf)
surfaceScalarField uS
(
    IOobject
    (
        "uS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    Uf & mesh.Sf()
);
uS.oldTime();
uS.write();

// The velocity in the dual perpendicular direction (the d direction)
surfaceScalarField v
(
    IOobject
    (
        "v", runTime.timeName(), dualMesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh.dualMap(Uf) & dualMesh.jdir()
);
v.oldTime();
v.write();

const volScalarField f
(
    IOobject("f", runTime.constant(), dualMesh),
    2*(dualMesh.Omega() & dualMesh.rHat())
);

volScalarField pv
(
    IOobject
    (
        "pv", runTime.timeName(), dualMesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    f + TRiSK3d::curl(v)
);
if (linear) pv = f;
pv.write();
surfaceScalarField pvf("pvf", fvc::interpolate(pv));

surfaceScalarField pvFlux("pvFlux", pvf*TRiSK3d::perp(uS));
if (conserveEnergy)
{
    pvFlux = 0.5*(pvFlux + TRiSK3d::perp(dualMesh.dualMap(pvf)*uS));
}

pvFlux.oldTime();

volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Phi", dimensionSet(0,2,-2,0,0), scalar(0)),
    "zeroGradient"
);

surfaceScalarField duSdt
(
    IOobject("duSdt", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("duSdt", dimensionSet(0,3,-2,0,0), scalar(0))
);
duSdt.oldTime();

