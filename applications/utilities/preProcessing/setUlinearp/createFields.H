// Read in prognostic variables

Info<< "Reading Exner\n" << endl;
const volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

const surfaceScalarField pf("pf", pRef*pow(fvc::interpolate(Exner), 1./kappa));

const surfaceVectorField Uf_init
(
    IOobject("Uf_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);


surfaceVectorField Uf
(
    IOobject("Uf", runTime.timeName(), mesh),
    Uf_init
);

