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
            fvo[faceI] = ?;
            fvn[faceI] = ?;
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
    
    // Scale the face vols so that they sum to each cell volume
    scalarField Vnew(nCells(), scalar(0));
    forAll(fvo, faceI)
    {
        Vnew[owner()[faceI]] += fvo[faceI];
        Vnew[neighbour()[faceI]] += fvn[faceI];
    }
    forAll(fvo.boundaryField(), patchI)
    {
        const labelUList& pFaceCells = boundary()[patchI].faceCells();
        forAll(fvo.boundaryField()[patchI], fi)
        {
            Vnew[pFaceCells[fi]] += fvo.boundaryField()[patchI][fi];
        }
    }
    
    forAll(fvo, faceI)
    {
        fvo[faceI] *= V()[owner()[faceI]]/Vnew[owner()[faceI]];
        fvn[faceI] *= V()[neighbour()[faceI]]/Vnew[neighbour()[faceI]];
    }
    forAll(fvo.boundaryField(), patchI)
    {
        const labelUList& pFaceCells = boundary()[patchI].faceCells();
        forAll(fvo.boundaryField()[patchI], fi)
        {
            fvo.boundaryField()[patchI][fi] *= V()[pFaceCells[fi]]
                                              /Vnew[pFaceCells[fi]];
        }
    }

    faceVolPtr_ = new surfaceScalarField
    (
        "faceVol",
        faceVolOwn() + faceVolNei()
    );
}

