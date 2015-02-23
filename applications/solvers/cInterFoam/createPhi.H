    // phi with boundaryies to ensure correct flow over boundaries
    wordList phiBCs(mesh.boundaryMesh().size(), "calculated");
    forAll(mesh.boundaryMesh(), patchi)
    {
        const word btype = mesh.boundaryMesh()[patchi].type();
        const word utype = U.boundaryField()[patchi].type();
        if
        (
            btype == "wall" || btype == "symmetryPlane" || btype == "empty"
         || utype == "fixedValue" || utype == "slip"
        )
        {
            phiBCs[patchi] = "fixedValue";
        }
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U) & mesh.Sf(),
        phiBCs
    );

    forAll(phi.boundaryField(), patchi)
    {
        const word btype = mesh.boundaryMesh()[patchi].type();
        const word utype = U.boundaryField()[patchi].type();
        if
        (
            btype == "wall" || btype == "symmetryPlane" || btype == "empty"
         || utype == "slip"
        )
        {
            phi.boundaryField()[patchi] == 0;
        }
    }
