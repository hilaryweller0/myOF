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
const dimensionedScalar beta(earthProperties.lookup("beta"));
const dimensionedScalar H = linear ?
                          earthProperties.lookup("H") :
                          dimensionedScalar("H", dimLength, scalar(0));
