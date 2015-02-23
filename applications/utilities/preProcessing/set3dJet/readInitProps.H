Info<< "Reading initial properties\n" << endl;

IOdictionary initProps
(
    IOobject
    (
        "initProps",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//const scalar etaTrop = readScalar(initProps.lookup("etaTrop"));
const dimensionedScalar u0(initProps.lookup("u0"));
//const scalar eta0 = readScalar(initProps.lookup("eta0"));
//const dimensionedScalar lapseRate(initProps.lookup("lapseRate"));
//const dimensionedScalar DT(initProps.lookup("DT"));

