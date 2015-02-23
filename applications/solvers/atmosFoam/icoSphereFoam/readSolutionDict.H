const scalar offCentre
(
    mesh.schemesDict().lookupOrDefault<scalar>("offCentre", scalar(0.5))
);
const bool linear
(
    mesh.solutionDict().lookupOrDefault<bool>("linear", "false")
);
const Switch conserveEnergy
(
    mesh.solutionDict().lookup("conserveEnergy")
);

// Read PISO controls
const dictionary& piso = mesh.solutionDict().subDict("PISO");
const int nCorr = piso.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr = piso.lookupOrDefault<int>("nNonOrthCorr", 1);

