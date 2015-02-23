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
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("h0", dimLength, scalar(0))
);
Info << "max h0 = " << max(h0) << endl;

const surfaceScalarField gsnGradh0 = magg*fvc::snGrad(h0);

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

surfaceVectorField UfD
(
    IOobject("UfD", runTime.timeName(), dualMesh),
    TRiSK::primalToDualFaceMap(Uf)
);

surfaceScalarField un
(
    IOobject("un", runTime.timeName(), mesh),
    (Uf & mesh.Sf())/mesh.magSf()
);
un.oldTime();

surfaceScalarField unD
(
    IOobject("unD", runTime.timeName(), dualMesh),
    TRiSK::primalToDualFluxMap(TRiSK::perp(un))
);

surfaceScalarField hf
(
    IOobject("hf", runTime.timeName(), mesh),
    fvc::interpolate(h)
//    TRiSK::cellToFaceMap(h)
);
if (linear) hf == H;
surfaceScalarField hS = hf*mesh.magSf();

volScalarField hfc("hfc", TRiSK::faceToCellMap(hf));
volScalarField hfcDiff("hfcDiff", hfc - h);
Info << "sum(h) = " << fvc::domainIntegrate(h).value()<<endl;
Info 
     << "\nsum(hf) = " << TRiSK::domainIntegrate(hf).value()<<endl;
Info
     << "\nsum(hfc) = " << fvc::domainIntegrate(hfc).value()<<endl;
Info
     << "\nnormalised difference 1,3 = "
     << ((fvc::domainIntegrate(h) - fvc::domainIntegrate(hfc))/fvc::domainIntegrate(h)).value()<<endl;
Info
     << "\nnormalised difference 1,2 = "
     << ((fvc::domainIntegrate(h) - TRiSK::domainIntegrate(hf))/fvc::domainIntegrate(h)).value()<<endl;
Info
     << "\nnormalised difference 2,3 = "
     << ((TRiSK::domainIntegrate(hf) - fvc::domainIntegrate(hfc))/fvc::domainIntegrate(hfc)).value()
 << endl;


surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh),
    un*hf*mesh.magSf()
);
phi.oldTime();

const volScalarField f
(
    IOobject("f", runTime.constant(), dualMesh),
    2*(Omega & dualMesh.C())/mag(dualMesh.C())
);

const surfaceScalarField ff("ff", 2*(Omega & mesh.Cf())/mag(mesh.Cf()));
/*(
    IOobject
    (
        "ff",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);*/
volScalarField hv(IOobject("hv", runTime.timeName(), dualMesh), dualMesh, H);
if (!linear) hv = TRiSK::primalToDualCellMap(h);

volScalarField vorticityByH
(
    IOobject("vorticityByH", runTime.timeName(), dualMesh),
    TRiSK::curl(un)/hv
);

volScalarField pv(IOobject("pv", runTime.timeName(), dualMesh), f/H);
if (!linear) pv = vorticityByH + f/hv;
pv.write();

surfaceScalarField pvf = TRiSK::dualToPrimalFaceMap(linearInterpolate(pv),mesh);

if (interpolateVorticity)
{
    pvf = ff/hf;
    if (!linear)
    {
        pvf += TRiSK::dualToPrimalFaceMap(fvc::interpolate(vorticityByH), mesh);
    }
}
else
{
    pvf = TRiSK::dualToPrimalFaceMap(fvc::interpolate(pv), mesh);
}

surfaceScalarField pvFlux = trisk ? pvf*TRiSK::perp(hf*un)
                                  : pvf*TRiSK::linPerp(hf*un);
pvFlux.oldTime();

if (conserveEnergy)
{
    if (trisk) pvFlux = 0.5*(pvFlux + TRiSK::perp(pvf*hf*un));
    else       pvFlux = 0.5*(pvFlux + TRiSK::linPerp(pvf*hf*un));
}


volScalarField KE
(
    IOobject("KE", runTime.timeName(), mesh),
    alphaKE*TRiSK::ke(un) + (1-alphaKE)*TRiSK::vertexKE(un)
);
KE.oldTime();

//surfaceScalarField dudt("dudt", -pvFlux - magg*fvc::snGrad(h) - gsnGradh0);
//dudt.oldTime();
//if (!linear) dudt -= fvc::snGrad(TRiSK::ke(un));

//Info << "pvFlux = " << pvFlux
//     << "\nend pvFlux" << endl;

