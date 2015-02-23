void fvMeshWithDual::makeSphFaceVols() const
{
    if (faceVolOwnPtr_ || faceVolNeiPtr_ || faceVolPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makefaceVols")
            << "faceVol already exists" << abort(FatalError);
    }

    faceVolOwnPtr_ = new surfaceScalarField
    (
        IOobject("faceVolOwn", pointsInstance(), *this),
        *this,
        dimensionedScalar("V", dimVol, scalar(0))
    );
    faceVolNeiPtr_ = new surfaceScalarField
    (
        IOobject("faceVolNei", pointsInstance(), *this),
        *this,
        dimensionedScalar("V", dimVol, scalar(0))
    );

    surfaceScalarField& fvo = *faceVolOwnPtr_;
    surfaceScalarField& fvn = *faceVolNeiPtr_;

    // for vertical internal faces
    forAll(fvo, faceI)
    {
        const label ie = faceToPatchEdge()[faceI];
        if (ie != -1)
        {
            // Set of faces making up the owner face volume
            const labelList& f = faces()[faceI];
            pointField fCtrs(f.size()+1);
            pointField fAreas(f.size()+1);
            fCtrs[f.size()] = Cf()[faceI];
            fAreas[f.size()] = Sf()[faceI];
            label iprev = f.size()-1;
            for(label i = 0; i < f.size(); i++)
            {
                pointField newFacePoints(4);
                newFacePoints[0] = mag(points()[f[i]])
                                  *unitVector(C()[owner()[faceI]]);
                newFacePoints[1] = points()[f[i]];
                newFacePoints[2] = points()[f[iprev]];
                newFacePoints[3] = mag(points()[f[iprev]])
                                  *unitVector(C()[owner()[faceI]]);
                calcSphFaceCentreAndArea(newFacePoints, fCtrs[i], fAreas[i]);
                iprev = i;
            }
                        
            fvo[faceI] = calcSphCellVol(fCtrs, fAreas);
            
            // Set of faces making up the neighbour face volume
            iprev = f.size()-1;
            for(label i = 0; i < f.size(); i++)
            {
                pointField newFacePoints(4);
                newFacePoints[0] = mag(points()[f[i]])
                                  *unitVector(C()[neighbour()[faceI]]);
                newFacePoints[1] = points()[f[i]];
                newFacePoints[2] = points()[f[iprev]];
                newFacePoints[3] = mag(points()[f[iprev]])
                                  *unitVector(C()[neighbour()[faceI]]);
                calcSphFaceCentreAndArea(newFacePoints, fCtrs[i], fAreas[i]);
                iprev = i;
            }
            
            fvn[faceI] = calcSphCellVol(fCtrs, fAreas);
        }
        else // horizontal face (ie == -1)
        {
            const point& fCtr = Cf()[faceI];
            const vector& fAr = Sf()[faceI];
            const point& Co = C()[owner()[faceI]];
            const point& Cn = C()[neighbour()[faceI]];
            fvo[faceI] = fAr & (fCtr - Co);
            fvn[faceI] = fAr & (Cn - fCtr);
        }
    }

    forAll(fvo.boundaryField(), patchI)
    {
        forAll(fvo.boundaryField()[patchI], fi)
        {
            const label faceI = boundaryMesh()[patchI].start() + fi;
            const label ie = faceToPatchEdge()[faceI];
            if (ie != -1) // vertical face
            {
                fvo.boundaryField()[patchI][fi] = 0;
            }
            else // horizontal face
            {
                fvo.boundaryField()[patchI][fi] = faceAreas()[faceI]
                     & (faceCentres()[faceI] - C()[faceOwner()[faceI]]);
            }
        }
    }

    forAll(fvn.boundaryField(), patchI)
    {
        forAll(fvn.boundaryField()[patchI], fi)
        {
            fvn.boundaryField()[patchI][fi] = 0;
        }
    }

//    // Scale the face vols so that they sum to each cell volume
//    scalarField Vnew(nCells(), scalar(0));
//    forAll(fvo, faceI)
//    {
//        Vnew[owner()[faceI]] += fvo[faceI];
//        Vnew[neighbour()[faceI]] += fvn[faceI];
//    }
//    forAll(fvo.boundaryField(), patchI)
//    {
//        const labelUList& pFaceCells = boundary()[patchI].faceCells();
//        forAll(fvo.boundaryField()[patchI], fi)
//        {
//            Vnew[pFaceCells[fi]] += fvo.boundaryField()[patchI][fi];
//        }
//    }
//    
////    Info << "Face vol difference = " << (Vnew - V())/V() << endl;

//    forAll(fvo, faceI)
//    {
//        fvo[faceI] *= V()[owner()[faceI]]/Vnew[owner()[faceI]];
//        fvn[faceI] *= V()[neighbour()[faceI]]/Vnew[neighbour()[faceI]];
//    }
//    forAll(fvo.boundaryField(), patchI)
//    {
//        const labelUList& pFaceCells = boundary()[patchI].faceCells();
//        forAll(fvo.boundaryField()[patchI], fi)
//        {
//            fvo.boundaryField()[patchI][fi] *= V()[pFaceCells[fi]]
//                                              /Vnew[pFaceCells[fi]];
//        }
//    }

    faceVolPtr_ = new surfaceScalarField
    (
        "faceVol",
        faceVolOwn() + faceVolNei()
    );
}

