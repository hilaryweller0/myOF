Info<< "Reading field h" << endl;
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

//Info<< "Reading field h0 if present" << endl;
//const volScalarField h0
//(
//    IOobject
//    (
//        "h0",
//        runTime.constant(),
//        mesh,
//        IOobject::READ_IF_PRESENT,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("h0", dimLength, scalar(0))
//);
//Info << "max h0 = " << max(h0) << endl;

//const surfaceScalarField gGradh0 = magg*mesh.magSf()
//                                  *mesh.signedDualMap(fvc::snGrad(h0));

Info << "Creating hf" << endl;
surfaceScalarField hf
(
    IOobject("hf", runTime.timeName(), mesh),
    mesh,
    H
);
if (!linear) hf = fvc::interpolate(h);

// The total velocity vector on the faces of the primal mesh
Info << "Reading Uf" << endl;
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

// The cell centre velocity
Info << "Reading U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Hops::cellReconstruct(Uf & mesh.Sf())
);
U.write();

// The flux in the dual perpendicular direction (the d direction)
Info << "Creating hvS" << endl;
surfaceScalarField hvS
(
    IOobject
    (
        "hvS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (Uf & mesh.ddir())*mesh.magSf()*hf
);
hvS.oldTime();
hvS.write();

Info<< "Creating Coriolis Acceleration" << endl;
const surfaceVectorField F
(
    "F",
    beta*mesh.Cf().component(vector::Y)*vector(0.,0.,1.)
);

//phi with boundaries to ensure correct flow over boundaries
wordList phiBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = h.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux" || htype == "zeroGradient"
    )
    {
        phiBCs[patchi] = "fixedValue";
    }
}

// The flux. Ie the velocity in the primal normal direction (hf*Sf)
Info << "Creating phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    (Uf & mesh.Sf())*hf,
    phiBCs
);
phi.oldTime();
phi.write();

forAll(phi.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = h.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux"
    )
    {
        phi.boundaryField()[patchi] == 0;
    }
}

Info << "Creating dvSdt" << endl;
surfaceScalarField dvSdt
(
    IOobject("dvSdt", runTime.timeName(), mesh),
    -mesh.magSf()*((mesh.ddir() & ((F ^ Uf)*hf)) + magg*hf*fvc::snGrad(h)),
    phi.boundaryField().types()
);
if (!linear)
{
    dvSdt -= mesh.magSf()*(mesh.ddir() & fvc::interpolate(fvc::div(phi,U)));
}
//dvSdt.write();
dvSdt.oldTime();

// reconstruct and U Uf from hvS to check code
surfaceVectorField Ufr
(
    IOobject("Ufr", runTime.timeName(), mesh),
    //Hops::faceReconstruct(Hops::ddirToFlux(hvS)/hf)
    Hops::faceReconstruct(hvS/hf, mesh.ddir())
);
Ufr.write();

volVectorField Ur
(
    IOobject("Ur", runTime.timeName(), mesh),
    //Hops::cellReconstruct(Hops::ddirToFlux(hvS)/hf)
    Hops::cellReconstruct(hvS/hf, mesh.ddir())
);
Ur.write();

Info << exit(FatalError);


