const Switch conserveEnergy
(
    mesh.solutionDict().lookup("conserveEnergy")
);

const Switch linear
(
    mesh.solutionDict().lookup("linear")
);

const Switch interpolateVorticity
(
    mesh.solutionDict().lookupOrDefault<Switch>("interpolateVorticity",false)
);

const Switch useTriskPerp
(
    mesh.solutionDict().lookupOrDefault<Switch>("useTriskPerp", true)
);

const Switch errorDiags
(
    mesh.solutionDict().lookupOrDefault<Switch>("errorDiags", true)
);

// Read PISO controls
const dictionary& piso = mesh.solutionDict().subDict("PISO");
const int nCorr = piso.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr = piso.lookupOrDefault<int>("nNonOrthCorr", 1);

const scalar offCentre
(
    mesh.schemesDict().lookupOrDefault<scalar>("offCentre", scalar(0.5))
);

// Read diffusion coefficients
const dimensionedScalar mu = mesh.solutionDict().lookupOrDefault<dimensionedScalar>
(
    "mu", dimensionedScalar("mu", dimensionSet(0,2,-1,0,0), scalar(0))
);
const Switch diffusion = (mag(mu.value()) > SMALL);
