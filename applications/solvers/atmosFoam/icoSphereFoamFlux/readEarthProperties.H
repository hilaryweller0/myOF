Info << "\nReading earthProperties" << endl;

IOdictionary earthProperties
(
    IOobject
    (
        "earthProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//const dimensionedScalar& magg = mesh.magg();
const dimensionedVector Omega = mesh.Omega();

