const scalar alphaKE
(
    mesh.solutionDict().lookupOrDefault<scalar>("alphaKE", scalar(1))
);
const scalar alphah
(
    mesh.solutionDict().lookupOrDefault<scalar>("alphah", scalar(1))
);

const Switch linear
(
    mesh.solutionDict().lookupOrDefault<Switch>("linear",false)
);

const Switch interpolateVorticity
(
    mesh.solutionDict().lookupOrDefault<Switch>("interpolateVorticity",false)
);

const Switch conserveEnergy
(
    mesh.solutionDict().lookupOrDefault<Switch>("conserveEnergy", true)
);

const Switch doubleAverage
(
    mesh.solutionDict().lookupOrDefault<Switch>("doubleAverage", false)
);

const Switch trisk
(
    mesh.solutionDict().lookupOrDefault<Switch>("trisk", true)
);

//const Switch newScheme
//(
//    mesh.solutionDict().lookupOrDefault<Switch>("newScheme", false)
//);

//const Switch colocateU
//(
//    mesh.solutionDict().lookupOrDefault<Switch>("colocateU", false)
//);

//const Switch fluxForm
//(
//    mesh.solutionDict().lookupOrDefault<Switch>("fluxForm", false)
//);

//const Switch replaceNormal
//(
//    mesh.solutionDict().lookupOrDefault<Switch>("replaceNormal", false)
//);

// Read PISO controls
const dictionary& piso = mesh.solutionDict().subDict("PISO");
const int nCorr = piso.lookupOrDefault<int>("nCorrectors", 1);

const scalar offCentre
(
    mesh.schemesDict().lookupOrDefault<scalar>("offCentre", scalar(0.5))
);


