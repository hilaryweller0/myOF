pointIOField newCellCentres
(
    IOobject("cellCentres", runTime.timeName(), "polyMesh", newMesh),
    diamondCentres
);

newCellCentres.write();

