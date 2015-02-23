const scalar Vtot = sum(mesh.V()).value();

// Sf in the k direction
const surfaceVectorField Sfhat = mesh.Sf()/mesh.magSf();
const surfaceScalarField kdir = mag(Sfhat.component(vector::Z));
const surfaceScalarField kSf = mesh.Sf().component(vector::Z);

const surfaceScalarField impDir = hevi ?
      surfaceScalarField("impDir", kdir) :
      surfaceScalarField("impDir", mag(Sfhat));

const surfaceScalarField expDir("expDir", -impDir + 1);

// Accoustic Courant number
{
    surfaceScalarField CohCoeff
    (
        "CohCoeff",
        runTime.deltaT()*mesh.deltaCoeffs()
          *mag(Sfhat - Sfhat.component(vector::Z)*vector(0,0,1))
    );
    surfaceScalarField CovCoeff
    (
        "CovCoeff",
        runTime.deltaT()*mesh.deltaCoeffs()*mag(Sfhat.component(vector::Z))
    );

    {
        Info << "Maximum horizontal accoustic Courant number = "
             << cs.value()*max(CohCoeff.internalField()) << endl;
        Info << "Maximum vertical accoustic Courant number = "
             << cs.value()*max(CovCoeff.internalField()) << endl;
    }
    Info << "dt*N = " << dt*N << endl;
}

Info << "Reading the stream function and taking the curl" << endl;
const surfaceScalarField streamFunc0
(
    IOobject("streamFunc0", runTime.constant(), mesh,IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("", dimensionSet(0,2,-2,0,0), scalar(0))
);
surfaceScalarField streamFunc
(
    IOobject
    (
        "streamFunc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    0*streamFunc0
);
surfaceScalarField curlSf
(
    IOobject("curlSf", runTime.constant(), mesh),
    fvc::interpolate
    (
        TRiSK3d::curl(streamFunc*vector(0,1,0))
    ) & mesh.Sf()
);
curlSf.write();

Info << "Reading u0" << endl;
surfaceScalarField u0
(
    IOobject("u0", runTime.constant(), mesh, IOobject::MUST_READ),
    mesh
);

// Read in prognostic variables

Info<< "Reading field un" << endl;
surfaceScalarField un
(
    IOobject
    (
        "un",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
if (runTime.value() == 0 && !linear) un = u0;
un.oldTime();

Info<< "Reading field P" << endl;
volScalarField P
(
    IOobject
    (
        "P",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
P.oldTime();

volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading/creating field bkSf" << endl;
surfaceScalarField bkSf
(
    IOobject
    (
        "bkSf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(b)*kSf
);
bkSf.oldTime();
volVectorField bk = fvc::reconstruct(bkSf);

// Creating diagnostic variables
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(un)
);
U.write();

surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    TRiSK3d::reconstructVec(un/mesh.magSf())
);
Uf.write();

const surfaceVectorField Ufbar
(
    IOobject("Ufbar", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
    0*Uf
);

const volVectorField Ubar
(
    IOobject("Ubar", runTime.constant(), mesh, IOobject::READ_IF_PRESENT),
    0*U
);

const volVectorField& Utot = linear? Ubar : U;
const surfaceVectorField& Uftot = linear? Ufbar : Uf;

const surfaceVectorField gradUbar = fvc::snGrad(Ubar);

const surfaceScalarField delta2 = 1./sqr(mesh.deltaCoeffs());

volScalarField divu
(
    IOobject("divu", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::div(un)
);
divu.write();

volScalarField divuNu
(
    IOobject("divuNu", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE),
    fvc::div(un)
);

