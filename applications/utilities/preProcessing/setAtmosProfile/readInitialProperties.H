Info << "\nReading initialProperties" << endl;

IOdictionary initialProperties
(
    IOobject
    (
        "initialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar BruntV(initialProperties.lookup("BruntVaisallaFreq"));

