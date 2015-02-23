// Create the mesh displacement potential, Phi
volScalarField Phi
(
    IOobject("Phi",runTime.timeName(),mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);
Phi.oldTime();

//volVectorField rHat = mesh.C()/mag(mesh.C());

// mimetic gradPhi
surfaceVectorField gradPhif
(
    IOobject("gradPhif", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    TRiSK::divReconstruct(fvc::snGrad(Phi)*mesh.magSf())
);
volVectorField gradPhi
(
    IOobject("gradPhi", runTime.timeName(), mesh),
    TRiSK::faceToCellMap(gradPhif)
);

//// non mimetic gradPhi
//volVectorField gradPhi
//(
//    IOobject("gradPhi", runTime.timeName(), mesh),
//    fvc::reconstruct(fvc::snGrad(Phi)*mesh.magSf())
//);

//surfaceVectorField gradPhif
//(
//    IOobject("gradPhif", runTime.timeName(), mesh, 
//             IOobject::NO_READ, IOobject::AUTO_WRITE),
//    fvc::interpolate(gradPhi)
//);

//gradPhif += (fvc::snGrad(Phi) - (gradPhif & mesh.Sf())/mesh.magSf())
//            *mesh.Sf()/mesh.magSf();

//// mimetic gradPhi
//surfaceVectorField gradPhif
//(
//    IOobject("gradPhif", runTime.timeName(), mesh, 
//             IOobject::NO_READ, IOobject::AUTO_WRITE),
//    TRiSK::divReconstruct(fvc::snGrad(Phi)*mesh.magSf())
//);

//volVectorField gradPhi
//(
//    IOobject("gradPhi", runTime.timeName(), mesh),
//    TRiSK::faceToCellMap(gradPhif)
//);

// determinant of the Hessian of Phi
volScalarField detHess
(
    IOobject("detHess", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("", dimless, scalar(1))
);
volTensorField Hessian("Hessian", fvc::grad(gradPhif));
volScalarField del2Phi
(
    IOobject("del2Phi", runTime.timeName(), mesh, 
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::laplacian(Phi)
);

// the volume ratio between old and new meshes
volScalarField volRatio
(
    IOobject("volRatio", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    rMesh, dimless, wordList(2,"zeroGradient") //"empty") //
);
volRatio.internalField() = 1;
volRatio.correctBoundaryConditions();
volRatio.write();

const volScalarField monitor
(
    IOobject("monitor", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar("monitor", dimensionSet(0,-2,0,0,0), scalar(1))
);

volScalarField equidist
(
    IOobject("equidist", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    monitor*detHess
);

volScalarField monitorNew("monitorNew", monitor);

volScalarField monitorR
(
    IOobject("monitorR", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    monitorFunc().map(rMesh, monitor)
);
monitorNew.internalField() = monitorR.internalField();

// New Voronoi generating points
pointIOField Vpoints
(
    IOobject("Vpoints", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh.C()
);

// Equidistribution constant
dimensionedScalar equiDistMean = fvc::domainIntegrate(detHess)
                                /fvc::domainIntegrate(1./monitorNew);

//// Source term to make the Poisson equation into the PMA equation
//volScalarField source
//(
//    IOobject("source", runTime.timeName(), mesh,
//             IOobject::NO_READ, IOobject::AUTO_WRITE),
//    equiDistMean/monitorNew + del2Phi - detHess
//);
//source.oldTime();

//// Scale the source term by this ratio
//volScalarField sourceRatio
//(
//    IOobject("sourceRatio", runTime.timeName(), mesh,
//             IOobject::NO_READ, IOobject::AUTO_WRITE),
//    source/source
//);
//sourceRatio.oldTime();
////scalar spread = GREAT;
////scalar spreadOld = GREAT;


// Stabilisation
volScalarField boostLaplacian
(
    IOobject("boostLaplacian", runTime.timeName(), rMesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    detHess/(1+del2Phi)
);

