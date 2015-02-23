const Switch conserveEnergy
(
    mesh.solutionDict().lookup("conserveEnergy")
);

const Switch linear
(
    mesh.solutionDict().lookup("linear")
);

const Switch incompressible
(
    mesh.solutionDict().lookupOrDefault<Switch>("incompressible", false)
);

// Read PISO controls
const dictionary& piso = mesh.solutionDict().subDict("PISO");
const int nCorr = piso.lookupOrDefault<int>("nCorrectors", 1);

