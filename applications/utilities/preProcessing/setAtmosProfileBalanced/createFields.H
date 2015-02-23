// Read in prognostic variables

Info<< "Creating Exner_init\n" << endl;
volScalarField Exner_init
(
    IOobject("Exner_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Setting Exner\n" << endl;
volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
    Exner_init
);

Info<< "Creating theta_init\n" << endl;
volScalarField theta_init
(
    IOobject("theta_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Setting theta\n" << endl;
volScalarField theta
(
    IOobject("theta", runTime.timeName(), mesh, IOobject::NO_READ),
    theta_init
);

surfaceScalarField thetaf("thetaf", fvc::interpolate(theta));

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
);
surfaceScalarField rhof = fvc::interpolate(rho);

surfaceScalarField phi
(
    IOobject("phi", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("phi", dimensionSet(1,0,-2,0,0), scalar(0))
);

