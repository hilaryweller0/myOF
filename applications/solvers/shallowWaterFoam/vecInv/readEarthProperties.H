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

const dimensionedScalar& magg = mesh.magg();
const dimensionedVector Omega = mesh.Omega();
const dimensionedScalar H = linear ? earthProperties.lookup("H") :
      dimensionedScalar("H", dimLength, scalar(1));

