Info<< "Reading field P\n" << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

//Boundary conditions to ensure correct flow over boundaries
wordList uBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = P.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux" || htype == "zeroGradient"
    )
    {
        uBCs[patchi] = "fixedValue";
    }
}

// The velocity flux in the face normal direction
surfaceScalarField uS
(
    IOobject
    (
        "uS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT
    ),
    Uf & mesh.Sf(),
    uBCs
);

forAll(uS.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = P.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux"
    )
    {
        uS.boundaryField()[patchi] == 0;
    }
}

// The velocity in the cell centre to cell centre direction *S
surfaceScalarField vS
(
    IOobject
    (
        "vS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT
    ),
    (Uf & mesh.ddir())*mesh.magSf()
);
vS.oldTime();

//surfaceScalarField vSd("vSd", -mesh.dualFluxMap(vS));

const volScalarField f
(
    IOobject("f", runTime.constant(), dualMesh),
    2*(Omega & dualMesh.rHat())
);
//const surfaceScalarField ff
//(
//    IOobject("ff", runTime.constant(), mesh),
//    2*(Omega & mesh.rHatf())
//);
//f.write();
//const volScalarField curlOmegax
//(
//    IOobject("curlOmegax", runTime.constant(), dualMesh),
//    TRiSK::curl((Omega ^ dualMesh.Cf()) & mesh.dualFluxMap(mesh.Sf()))
//);
//curlOmegax.write();
//FatalErrorIn("") << exit(FatalError);

volScalarField pv
(
    IOobject("pv", runTime.timeName(), dualMesh,
              IOobject::NO_READ, IOobject::AUTO_WRITE),
    f
);
if (!linear) pv += TRiSK::curl(mesh.dualFluxMap(vS));
surfaceScalarField pvzf = dualMesh.dualMap(fvc::interpolate(pv));

surfaceVectorField pvf
(
    IOobject
    (
        "pvf", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    2*Omega + TRiSK::curl3d(vS)
);
if (linear) pvf = 2*Omega;

surfaceVectorField pvh("pvh", pvf - (pvf & mesh.rHatf())*mesh.rHatf());

surfaceScalarField pvFlux
(
    "pvFlux",
    (pvf ^ Uf) & mesh.Sf()
);

volScalarField KE
(
    IOobject
    (
        "KE", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    TRiSK::KE(uS, vS)
);
KE.write();

surfaceVectorField djperp("djperp", unitVector(mesh.ddir() ^ mesh.jdir()));
//surfaceScalarField kk("kk", djperp & mesh.rHatf());
//kk.write();
#include "velocityVorticity.H"

Info << "Creating dvSdt" << endl;
surfaceScalarField dvSdt
(
    IOobject("dvSdt", runTime.timeName(), mesh),
    -pvFlux - mesh.magSf()*fvc::snGrad(P)
);
if (!linear)
{
    dvSdt -= mesh.magSf()*fvc::snGrad(KE);
}
dvSdt.oldTime();


