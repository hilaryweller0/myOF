// Adjust the Hessian to account for the anisotropic cells

const surfaceScalarField deltas = 1/rMesh.deltaCoeffs();

forAll(detHess, cellI)
{
    // Find minDx and maxDx for this cell
    const labelList& cellFaces = mesh.cells()[cellI];
    scalar minDx = GREAT;
    scalar maxDx = 0;
    for(label fi = 0; fi < cellFaces.size(); fi++)
    {
        label faceI = cellFaces[fi];
        if (faceI < mesh.nInternalFaces())
        {
            scalar dx = deltas[faceI];
            if (dx <= minDx) minDx = dx;
            if (dx >= maxDx) maxDx = dx;
        }
    }

    // set the volume to maxDx^2 or minDx^2 depending on detHess
    if (detHess[cellI] < volRatio[cellI])
    {
        detHess[cellI] = sqr(maxDx)/mesh.V()[cellI];
    }
    else if (detHess[cellI] > volRatio[cellI])
    {
        detHess[cellI] = sqr(minDx)/mesh.V()[cellI];
    }
}

