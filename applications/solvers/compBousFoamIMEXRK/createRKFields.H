Info << "Reading Butcher Tableaus\n" << endl;
const butcherTableau exBt
(
    mesh.schemesDict().subDict("ddtSchemes").lookup
    (
        "explicitButcherTableau"
    )
);
const butcherTableau imBt
(
    mesh.schemesDict().subDict("ddtSchemes").lookup
    (
        "implicitButcherTableau"
    )
);

// Fields to be stored at each Runge-Kutta sup-step

RKfield<scalar, fvsPatchField, surfaceMesh> sdBdt
(
    exBt.nSteps(),
    0*bkSf/dt
);

RKfield<scalar, fvsPatchField, surfaceMesh> fdBdt
(
    imBt.nSteps(),
    0*bkSf/dt
);

RKfield<scalar, fvsPatchField, surfaceMesh> sdudt
(
    exBt.nSteps(),
    0*un/dt
);

RKfield<scalar, fvsPatchField, surfaceMesh> fdudt
(
    imBt.nSteps(),
    0*un/dt
);

RKfield<scalar, fvPatchField, volMesh> sdPdt
(
    exBt.nSteps(),
    0*P/dt
);

RKfield<scalar, fvPatchField, volMesh> fdPdt
(
    imBt.nSteps(),
    0*P/dt
);

// scale factor for treatment of semi-implicit gravity waves
RKfield<scalar, fvsPatchField, surfaceMesh> G
(
    imBt.nSteps(),
    surfaceScalarField("G", mag(Sfhat))
);

if (SIgwaves)
{
    for(label iRK = 0; iRK < imBt.nSteps(); iRK++)
    {
        G[iRK] = 1./(1.+sqr(dt*N*aii)*kdir);
    }
}

