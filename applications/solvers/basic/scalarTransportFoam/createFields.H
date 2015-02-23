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
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Calculate phi and make it divergence free
surfaceScalarField phi("phi", Uf & mesh.Sf());

volScalarField p
(
    IOobject("p", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("p", dimensionSet(0,2,-1,0,0), scalar(0)),
    "zeroGradient"
);
fvScalarMatrix pEqn
(
    fvm::laplacian(p) + fvc::div(phi)
);
pEqn.setReference(0, scalar(0));
pEqn.solve();

phi += pEqn.flux();

Uf += (phi*mesh.Sf() - (Uf & mesh.Sf())*mesh.Sf())/sqr(mesh.magSf());
Uf.write();

