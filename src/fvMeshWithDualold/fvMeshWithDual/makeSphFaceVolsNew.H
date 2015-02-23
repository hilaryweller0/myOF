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

    // for all vertial faces
    forAll(fvo, faceI)
    {const label ie = faceToPatchEdge()[faceI]; if (ie != -1)
    {
        // Identify the spherical faces above and below the owner and nei
        const label own = owner()[faceI];
        const label nei = neighbour()[faceI];
        label faceUpOwn = -1;
        label faceDownOwn = -1;
        label faceUpNei = -1;
        label faceDownNei = -1;
        const labelList& facesOwn = cells()[own];
        const labelList& facesNei = cells()[nei];
        // for the owner cell
        for(label i = 0; i < facesOwn.size() && faceDownOwn == -1; i++)
        {
            label faceII = facesOwn[i];
            if (dualFaceMap()[faceII] == -1)
            {
                if (faceUpOwn == -1) faceUpOwn = faceII;
                else faceDownOwn = faceII;
            }
        }
        if (faceUpOwn == -1 || faceDownOwn == -1)
        {
            FatalErrorIn("fvMeshWithDual::makeSphFaceVols")
                << " cannot find faces above and below cell " << own
                << exit(FatalError);
        }
        // Swap faceUpOwn and faceDownOwn if necessary
        if (magSqr(Cf()[faceDownOwn]) > magSqr(Cf()[faceUpOwn]))
        {
            label faceTmp = faceUpOwn;
            faceUpOwn = faceDownOwn;
            faceDownOwn = faceTmp;
        }

        // for the neighbour cell
        for(label i = 0; i < facesNei.size() && faceDownNei == -1; i++)
        {
            label faceII = facesNei[i];
            if (dualFaceMap()[faceII] == -1)
            {
                if (faceUpNei == -1) faceUpNei = faceII;
                else faceDownNei = faceII;
            }
        }
        if (faceUpNei == -1 || faceDownNei == -1)
        {
            FatalErrorIn("fvMeshWithDual::makeSphFaceVols")
                << " cannot find faces above and below cell " << nei
                << exit(FatalError);
        }
        // Swap faceUpNei and faceDownNei if necessary
        if (magSqr(Cf()[faceDownNei]) > magSqr(Cf()[faceUpNei]))
        {
            label faceTmp = faceUpNei;
            faceUpNei = faceDownNei;
            faceDownNei = faceTmp;
        }
        
        // Find the lower triangle
        // circulate around the face until we are going downwards
        bool goingDown = false;
        const labelList& face = faces()[faceI];
        label iThis = 0;
        label iPrev = face.size()-1;
        for(iThis = 0; iThis < face.size() && !goingDown; iThis++)
        {
            const point& pt0 = points()[face[iPrev]];
            const point& pt1 = points()[face[iThis]];
            if
            (
                mag(1 - (pt0 & pt1)/(mag(pt0)*mag(pt1))) < 0.1
             && magSqr(pt1) < magSqr(pt0)
            )
            {
                goingDown = true;
                iPrev = iThis+1;
                if (iPrev >= face.size()) iPrev = 0;
            }
            else iPrev = iThis;
        }
        
        // Lower triangle for the owner
        point p0 = points()[face[iThis]];
        point p1 = points()[face[iPrev]];
        point p2 = Cf()[faceDownOwn];
        
        // Find face centre of triangle below the owner cell
        point CfdownOwn = (p0 + p1 + p2)/3.;
        
        // the face centre of the triangle below the neighbour cell
        p2 = Cf()[faceDownNei];
        point CfdownNei = (p0 + p1 + p2)/3.;
        
        // Find the upper triangle
        // circulate around the face until we are going upwards
        bool goingUp = false;
        iThis = 0;
        iPrev = face.size()-1;
        for(iThis = 0; iThis < face.size() && !goingUp; iThis++)
        {
            const point& pt0 = points()[face[iPrev]];
            const point& pt1 = points()[face[iThis]];
            if
            (
                mag(1 - (pt0 & pt1)/(mag(pt0)*mag(pt1))) < 0.1
             && magSqr(pt1) > magSqr(pt0)
            )
            {
                goingUp = true;
                iPrev = iThis+1;
                if (iPrev >= face.size()) iPrev = 0;
            }
            else iPrev = iThis;
        }
        
        // Upper triangle of the owner
        p0 = points()[face[iThis]];
        p1 = points()[face[iPrev]];
        p2 = Cf()[faceUpOwn];
        
        // Find face centre of triangle above of the owner
        point CfupOwn = (p0 + p1 + p2)/3.;
        
        // upper triangle of the neighbour
        p2 = Cf()[faceUpNei];
        
        // Face centre of the triangle above the neighbour
        point CfupNei = (p0 + p1 + p2)/3.;
        
        fvo[faceI] = edgeAreaOwn()[ie]*mag(mag(CfupOwn)-mag(CfdownOwn))
           * sqr(0.5*(mag(CfupOwn) + mag(CfdownOwn)));
//         * 0.5*(mag(faceAreas()[faceUpOwn]) + mag(faceAreas()[faceDownOwn]));
        fvn[faceI] = edgeAreaNei()[ie]*mag(mag(CfupNei)-mag(CfdownNei))
           * sqr(0.5*(mag(CfupNei) + mag(CfdownNei)));
//         * 0.5*(mag(faceAreas()[faceUpNei]) + mag(faceAreas()[faceDownNei]));
    }}

    label celli = 0;
    Info << "cell " << celli << " has volume " << V()[celli] << " and faces\n";
    const labelList& fs = cells()[celli];
    Info << fs << "\n with volumes\n";
    scalar sum = 0;
    forAll(fs, facei)
    {
        label faceI = fs[facei];
        if (faceI < nInternalFaces())
        {
            if (owner()[faceI] == celli)
            {
                sum += faceVolOwn()[faceI];
                Info << faceVolOwn()[faceI] << ' ';
            }
            else if (neighbour()[faceI] == celli)
            {
                sum += faceVolNei()[faceI];
                Info << faceVolNei()[faceI] << ' ';
            }
        }
    }
    Info << nl << " which sum to " << sum << endl;

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

