// Set the HessianVolumeRatio based on the aspect ratio of each cell
const surfaceScalarField deltas = 1/rMesh.deltaCoeffs();
forAll(HessianVolumeRatio, cellI)
{
//    if (equidist[cellI] < equiDistMean.value()) HessianVolumeRatio[cellI] = 1;
//    else
//    {
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
                if(dx >= maxDx) maxDx = dx;
            }
        }
        aspectRatio[cellI] = max(maxDx/minDx/aspectRatio0[cellI], scalar(1));
        //aspectRatio[cellI] = max(maxDx/minDx, maxDx/minDx/aspectRatio0[cellI]);
        //aspectRatio[cellI] = maxDx/minDx;
        //aspectRatio[cellI] = max(maxDx/minDx, aspectRatio0[cellI]);
        //aspectRatio[cellI] = 1 + max(maxDx/minDx - aspectRatio0[cellI], scalar(0));
        HessianVolumeRatio[cellI] = Foam::pow(aspectRatio[cellI], scalar(-0.04));
//    }
}

// Smooth HessianVolumeRatio
//{
//    HessianVolumeRatio += 1e-3*fvc::laplacian
//    (
//        1/sqr(mesh.deltaCoeffs()),
//        HessianVolumeRatio
//    );

//    volScalarField HessianVolumeRatioSmooth = HessianVolumeRatio;
//    forAll(HessianVolumeRatio, cellI)
//    {
//        const labelList& CCs = mesh.cellCells()[cellI];
//        forAll(CCs, ci)
//        {
//            label cellII = CCs[ci];
//            HessianVolumeRatioSmooth[cellI] += HessianVolumeRatio[cellII];
//        }
//        HessianVolumeRatioSmooth[cellI] /= CCs.size()+1;
//    }
//    HessianVolumeRatio = HessianVolumeRatioSmooth;
//}
