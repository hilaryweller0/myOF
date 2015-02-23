// Read in prognostic variables

Info<< "Reading T\n" << endl;
const volScalarField T
(
    IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ),
    mesh
);

Info<< "Reading Exner_init\n" << endl;
const volScalarField Exner_init
(
    IOobject("Exner_init", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

//wordList ExnerBCs(Exner_init.boundaryField().types());
label groundBC = -1;
label topBC = -1;
forAll(mesh.boundaryMesh(), patchi)
{
    if (mesh.boundaryMesh()[patchi].name() == "top")
    {
//        ExnerBCs[patchi] = "fixedValue";
        topBC = patchi;
    }
    else if (mesh.boundaryMesh()[patchi].name() == "ground")
    {
        groundBC = patchi;
    }
}
if (groundBC == -1 || topBC == -1)
{
    FatalErrorIn("setExnerThetaFromT_H")
        << " no boundary called top or ground. The boundaries are called "
        << mesh.boundaryMesh().names()
        << exit(FatalError);
}

Info<< "Setting Exner\n" << endl;
volScalarField Exner
(
    IOobject("Exner", runTime.timeName(), mesh, IOobject::NO_READ),
    Exner_init
//    ExnerBCs
);

Info<< "Reading theta_init\n" << endl;
const volScalarField theta_init
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

//Info<< "Creating field rho\n" << endl;
//volScalarField rho
//(
//    IOobject("rho", runTime.timeName(), mesh),
//    pRef/(R*theta)*pow(Exner, (1-kappa)/kappa)
//);
//surfaceScalarField rhof = fvc::interpolate(rho);

surfaceScalarField U
(
    IOobject("U", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("U", dimensionSet(0,3,-2,0,0), scalar(0))
);

