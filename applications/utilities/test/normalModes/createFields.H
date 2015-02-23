Info<< "Creating field Uf\n" << endl;
surfaceVectorField Uf
(
    IOobject("Uf", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedVector("Uf", dimVelocity, vector(0,0,0))
);

Info<< "Creating field h\n" << endl;
volScalarField h
(
    IOobject("h", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("h", dimLength, scalar(0))
);

