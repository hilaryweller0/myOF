    dictionary MCSPPA = mesh.solutionDict().subDict("MCSPPA");

    Switch MCS(MCSPPA.lookup("MCS"));

    int nConvectionSteps = 0;

    if (MCS)
    {
        MCSPPA.lookup("nConvectionSteps") >> nConvectionSteps;
    }

    int nNonOrthCorr = 0;
    if (MCSPPA.found("nNonOrthogonalCorrectors"))
    {
        MCSPPA.lookup("nNonOrthogonalCorrectors") >> nNonOrthCorr;
    }


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, MCSPPA, pRefCell, pRefValue);
