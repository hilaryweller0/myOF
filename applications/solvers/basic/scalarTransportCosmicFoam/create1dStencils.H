// Read in a uniform mesh and create 1d stencils using the points from this mesh

Info << "Creating uniform mesh for time " << runTime.timeName() << endl;

fvMesh meshU
(
    Foam::IOobject
    (
        "uniformMesh",
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

// Save points from normal mesh
pointField meshPoints = mesh.points();

// Move the points of the normal meshes in order to create 1d stencils
mesh.movePoints(meshU.points());

// create stencils for 1d interpolations
volVectorField gradT = fvc::grad(fvc::interpolate(T));

// Move the mesh points back again
mesh.movePoints(meshPoints);

