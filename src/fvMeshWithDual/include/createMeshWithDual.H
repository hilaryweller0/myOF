Foam::Info
    << "Create mesh with dual for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::word regionName;

if (args.optionReadIfPresent("region", regionName))
{
    Foam::Info
        << "Create mesh " << regionName << " for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}
else
{
    regionName = Foam::fvMesh::defaultRegion;
    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}

Foam::fvMeshWithDual mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    true,
    fvMeshWithDual::SPHERICALDIST
);

