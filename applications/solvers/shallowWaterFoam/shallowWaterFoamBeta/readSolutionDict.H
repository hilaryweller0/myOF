// Read PISO controls
const dictionary& piso = mesh.solutionDict().subDict("PISO");
const int nCorr = piso.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr = piso.lookupOrDefault<int>("nNonOrthCorr", 1);

const scalar offCentre
(
    mesh.schemesDict().lookupOrDefault<scalar>("offCentre", scalar(0.5))
);
const bool linear
(
    mesh.schemesDict().lookupOrDefault<bool>("linear", "false")
);

