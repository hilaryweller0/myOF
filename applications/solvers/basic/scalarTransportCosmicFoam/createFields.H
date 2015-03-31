Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field Uf\n" << endl;

surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Velocity in each direction
volScalarField u("u", U.component(vector::X));
//volScalarField v("v", Uf.component(vector::Y));
volScalarField w("w", U.component(vector::Z));
surfaceScalarField uf("uf", Uf.component(vector::X));
surfaceScalarField wf("wf", Uf.component(vector::Z));

// and the flux for the Courant number
surfaceScalarField phi("phi", Uf & mesh.Sf());

