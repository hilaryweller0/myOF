Info<< "Reading field P\n" << endl;
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

// The total velocity vector on the faces of the primal mesh
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// The cell centre velocity
Info << "Reading U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Hops::cellReconstruct(Uf & mesh.Sf())
);
U.write();

//Boundary conditions to ensure correct flow over boundaries
wordList uBCs(mesh.boundaryMesh().size(), "calculated");
forAll(mesh.boundaryMesh(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = P.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux" || htype == "zeroGradient"
    )
    {
        uBCs[patchi] = "fixedValue";
    }
}

// The velocity flux in the face normal direction
surfaceScalarField uS
(
    IOobject
    (
        "uS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT
    ),
    Uf & mesh.Sf(),
    uBCs
);

forAll(uS.boundaryField(), patchi)
{
    const word btype = mesh.boundaryMesh()[patchi].type();
    const word htype = P.boundaryField()[patchi].type();
    if
    (
        btype == "wall" || btype == "symmetryPlane" || btype == "empty"
     || htype == "geostrophicZeroFlux"
    )
    {
        uS.boundaryField()[patchi] == 0;
    }
}

// The velocity in the cell centre to cell centre direction *S
surfaceScalarField vS
(
    IOobject
    (
        "vS", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT
    ),
    (Uf & mesh.ddir())*mesh.magSf()
);
vS.oldTime();

Info << "Creating dvSdt" << endl;
surfaceScalarField dvSdt
(
    IOobject("dvSdt", runTime.timeName(), mesh),
    -mesh.magSf()*((mesh.ddir() & (Omega ^ Uf)) + fvc::snGrad(P)),
    uS.boundaryField().types()
);
if (!linear)
{
    dvSdt -= mesh.magSf()*(mesh.ddir() & fvc::interpolate(fvc::div(uS,U)));
}
//dvSdt.write();
dvSdt.oldTime();


