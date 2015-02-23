Info<< "Creating field g.h at cell centres\n" << endl;
volScalarField gh("gh", g & mesh.C());

Info<< "Creating Exner_init\n" << endl;
volScalarField Exner_init
(
    IOobject("Exner_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Creating theta_init\n" << endl;
volScalarField theta_init
(
    IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

