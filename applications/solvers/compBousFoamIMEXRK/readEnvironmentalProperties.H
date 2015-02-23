Info << "\nReading environmentalProperties" << endl;

IOdictionary envProps
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar omega(envProps.lookup("omega"));

const dimensionedScalar N(envProps.lookup("N"));
const dimensionedScalar cs(envProps.lookup("speedSound"));

const dimensionedScalar N2 = sqr(N);
const dimensionedScalar cs2 = sqr(cs);

const dimensionedScalar K
(
    envProps.lookupOrDefault<dimensionedScalar>
    (
        "K",
        dimensionedScalar("", dimless, scalar(0))
    )
);

const bool hyperDiff = K.value() > VSMALL;

