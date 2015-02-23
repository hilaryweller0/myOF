/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meshWithDual.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //


tmp<surfaceScalarField> fvMeshWithDual::magDelta() const
{
    if (debug)
    {
        Info<< "void fvMeshWithDual::magDelta() : "
            << "calculating face deltas"
            << endl;
    }

    tmp<surfaceScalarField> tdelta
    (
        new surfaceScalarField
        (
            IOobject
            (
                "magDelta",
                pointsInstance(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimLength
        )
    );
    surfaceScalarField& delta = tdelta();

    const volVectorField& C = this->C();
    const labelUList& owner = this->owner();
    const labelUList& neighbour = this->neighbour();

    if (isSpherical())
    {
        forAll(owner, facei)
        {
            delta[facei] = sphDist(C[neighbour[facei]], C[owner[facei]]);
        }
    }
    else
    {
        forAll(owner, facei)
        {
            delta[facei] = mag(C[neighbour[facei]] - C[owner[facei]]);
        }
    }

    forAll(delta.boundaryField(), patchi)
    {
        delta.boundaryField()[patchi] = mag(boundary()[patchi].delta());
    }

    return tdelta;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshWithDual::readFaceToPatchEdge() const
{
    if (faceToPatchEdgePtr_) FatalErrorIn("fvMeshWithDual::readFaceToPatchEdge")
                         << "faceToPatchEdge already exists" << abort(FatalError);
    faceToPatchEdgePtr_ = new IOList<label>
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readFaceToPatchFace() const
{
    if (faceToPatchFacePtr_) FatalErrorIn("fvMeshWithDual::readFaceToPatchFace")
                         << "faceToPatchFace already exists" << abort(FatalError);
    faceToPatchFacePtr_ = new IOList<label>
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readPointToPatchPoint() const
{
    if(pointToPatchPointPtr_)FatalErrorIn("fvMeshWithDual::readPointToPatchPoint")
                        << "pointToPatchPoint already exists" << abort(FatalError);
    pointToPatchPointPtr_ = new IOList<label>
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readCellToPatchFace() const
{
    if (cellToPatchFacePtr_) FatalErrorIn("fvMeshWithDual::readCellToPatchFace")
                         << "cellToPatchFace already exists" << abort(FatalError);
    cellToPatchFacePtr_ = new IOList<label>
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}

void Foam::fvMeshWithDual::readDualFaceMap() const
{
    if (dualFaceMapPtr_) FatalErrorIn("fvMeshWithDual::readDualFaceMap")
                         << "dualFaceMap already exists" << abort(FatalError);
    dualFaceMapPtr_ = new IOList<label>
    (
        IOobject
        (
            "dualFaceMap",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    );
}


void Foam::fvMeshWithDual::makesignFaceMap() const
{
    if (signFaceMapPtr_) FatalErrorIn("fvMeshWithDual::makesignFaceMap")
                           << "signFaceMap already exists" << abort(FatalError);
    signFaceMapPtr_ = new List<bool>(dualFaceMap().size());
    List<bool>& sfm = *signFaceMapPtr_;

    forAll(dualFaceMap(), faceI)
    {
        if (dualFaceMap()[faceI] != -1)
        {
            label faced = dualFaceMap()[faceI];
            sfm[faceI] = sign
            (
                (faceAreas()[faceI] ^ dualMesh().faceAreas()[faced])
              & faceCentres()[faceI]
            ) == 1;
        }
    }
}

//#include "makeSphFaceVolsNew.H"
#include "makeSphFaceVolsOld.H"
//#include "makeSphFaceVols210.H"

void fvMeshWithDual::makeFaceDepth() const
{
    if (depthfPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeFaceDepth")
            << " depthfPtr_ already exists" << exit(FatalError);
    }
    depthfPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "depthf",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this
        ),
        *this,
        dimensionedScalar("depthf", dimLength, scalar(0))
    );
    surfaceScalarField& d = *depthfPtr_;
    
    forAll(d, faceI)
    {
        // Find the edge of the patch for this (vertical) face
        const label ie = faceToPatchEdge()[faceI];
        if (ie != -1)
        {
            // points at either end of the edge
            const point p0 = unitVector
            (
                bottomPatch_.localPoints()[bottomPatch_.edges()[ie][0]]
            );
            const point p1 = unitVector
            (
                bottomPatch_.localPoints()[bottomPatch_.edges()[ie][1]]
            );
            
            const scalar faceLength = mag(Cf()[faceI])*sphDist(p0,p1);
            d[faceI] = magSf()[faceI]/faceLength;
        }
    }
}

void fvMeshWithDual::makeSphAreaFracs() const
{
    if (fpeAreaFracPtr_ || facePointAreaFracPtr_ || facePointsPtr_
        || edgeAreaOwnPtr_ || edgeAreaNeiPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeSphAreaFracs")
            << "fpeAreaFracPtr_ || facePointAreaFracPtr_ || facePointsPtr_"
            << " already exist" << exit(FatalError);
    }
    
    fpeAreaFracPtr_= new List<List<FixedList<scalar,2> > >(bottomPatch_.size());
    facePointAreaFracPtr_ = new scalarListList(bottomPatch_.size());
    facePointsPtr_ = new labelListList(bottomPatch_.size());
    edgeAreaOwnPtr_ = new scalarList(bottomPatch_.nEdges());
    edgeAreaNeiPtr_ = new scalarList(bottomPatch_.nEdges());
    
    List<List<FixedList<scalar,2> > >& fpeAfrac = *fpeAreaFracPtr_;
    scalarListList& fpAfrac = *facePointAreaFracPtr_;
    labelListList& fPts = *facePointsPtr_;
    scalarList& edgeAo = *edgeAreaOwnPtr_;
    scalarList& edgeAn = *edgeAreaNeiPtr_;

    // Set the patch edgeCentres from the vertical face centres
    pointField Ce(bottomPatch_.nEdges());
    forAll(Cf(), faceI)
    {
        label ie = faceToPatchEdge()[faceI];
        if (ie != -1) Ce[ie] = intersections()[faceI];
    }
    Ce = unitVector(Ce);

    // Circulate around each point of each face and calculate the proportion
    // of the face area associated with each point and each edge
    forAll(bottomPatch_, faceI)
    {
        // A list of edges for this face (for consistency with the stencil)
        const labelList& faceEdges = bottomPatch_.faceEdges()[faceI];
        
        const face& f = bottomPatch_[faceI];
        const point Ci = unitVector(bottomPatch_.faceCentres()[faceI]);
        
        // The dual point in this face
        //!! Warning, this won't word if either mesh is renumbered!!?!?
        //!! Warning, this won't word if either mesh is renumbered!!?!?
        //const point Ci = unitVector(dualMesh().bottomPatch().points()[faceI]);
        //!! Warning, this won't word if either mesh is renumbered!!?!?
        //!! Warning, this won't word if either mesh is renumbered!!?!?
        
        // Initialise the area fraction arrays and the facePoints array
        fpeAfrac[faceI].setSize(f.size());
        fpAfrac[faceI].setSize(f.size());
        fpAfrac[faceI] = scalar(0);
        fPts[faceI].setSize(f.size());
        
        // total for scaling to make into fractions
        scalar Atot = 0;
        
        // Circulate around face and calculate the half areas for each point
        label prevIp = bottomPatch_.edges()[faceEdges[0]][0];
        if
        (
            prevIp == bottomPatch_.edges()[faceEdges[1]][0]
         || prevIp == bottomPatch_.edges()[faceEdges[1]][1]
        )
        {
            prevIp = bottomPatch_.edges()[faceEdges[0]][1];
        }
        label prevIe = faceEdges.size()-1;
        forAll(faceEdges, iee)
        {
            const label ie = faceEdges[iee];
            label ip = bottomPatch_.edges()[ie][1];
            if (ip == prevIp) ip = bottomPatch_.edges()[ie][0];

            // current and previous points
            const point thisPoint = unitVector(bottomPatch_.localPoints()[ip]);
            const point prevPoint = unitVector(bottomPatch_.localPoints()[prevIp]);
            
            const scalar Ap = sphTriSolidAngle(Ce[ie], Ci, prevPoint);
            const scalar At = sphTriSolidAngle(Ce[ie], Ci, thisPoint);
            //const scalar A = sphTriSolidAngle(Ci, prevPoint, thisPoint);
//            Info << ie << " " << iee << " points\n" << Ce[ie] << nl << Ci
//                 << nl << prevPoint << nl << thisPoint << "\nareas\n"
//                 << Ap << nl << At << nl << A << nl;

            fpeAfrac[faceI][iee][0]    = At;
            fpeAfrac[faceI][prevIe][1] = Ap;
            fpAfrac[faceI][iee] += At;
            fpAfrac[faceI][prevIe] += Ap;
            fPts[faceI][iee] = ip;
            Atot += At + Ap;

            if (faceI == bottomPatch_.edgeFaces()[ie][0])
            {
                edgeAo[ie] = At+Ap;
            }
            else
            {
                edgeAn[ie] = At+Ap;
            }
            
            prevIp = ip;
            prevIe = iee;
        }

        // scale to make into fractions
        forAll(faceEdges, iee)
        {
            fpeAfrac[faceI][iee][0] /= Atot;
            fpeAfrac[faceI][iee][1] /= Atot;
            fpAfrac[faceI][iee] /= Atot;
//            const label ie = faceEdges[iee];
//            if (faceI == bottomPatch_.edgeFaces()[ie][0])
//            {
//                edgeAo[ie]/= Atot;
//            }
//            else
//            {
//                edgeAn[ie]/= Atot;
//            }
        }
    }

}

void fvMeshWithDual::makeIntersections() const
{
    if (intersectionsPtr_) FatalErrorIn("fvMeshWithDual::makeIntersections")
                        << "intersections already exists" << abort(FatalError);

    intersectionsPtr_ = new surfaceVectorField("intersections", Cf());

    surfaceVectorField& ints = *intersectionsPtr_;

    forAll(ints, faceI)
    {
        label faced = dualFaceMap()[faceI];
        if (faced != -1)
        {
//            const point& Cpo = C()[owner()[faceI]];
//            const point& Cpn = C()[neighbour()[faceI]];
//            const point& Cdo = dualMesh().C()[dualMesh().faceOwner()[faced]];
//          const point& Cdn = dualMesh().C()[dualMesh().faceNeighbour()[faced]];
//            plane pFace(Cpo^Cpn);
//            plane dFace(Cdo^Cdn);

            // Plane of faceI
            plane pFace
            (
                points()[faces()[faceI][0]],
                points()[faces()[faceI][1]],
                points()[faces()[faceI][2]]
            );
            
            // Plane of faced
            plane dFace
            (
                dualMesh().points()[dualMesh().faces()[faced][0]],
                dualMesh().points()[dualMesh().faces()[faced][1]],
                dualMesh().points()[dualMesh().faces()[faced][2]]
            );
            
            // intersection point
            plane::ray r = pFace.planeIntersect(dFace);
            ints[faceI] = unitVector(r.dir())*sign(r.dir() & Cf()[faceI])
                          *mag(ints[faceI]);
        }
    }
}

void fvMeshWithDual::makeCellEdgeInfo() const
{
    if (cellEdgeVolsPtr_ || cellEdgeCellsPtr_)
    {
        FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
             << "cellEdgeVolsPtr or cellEdgeCellsPtr already exists"
             << abort(FatalError);
    }
    
    cellEdgeVolsPtr_ = new scalarListList(nCells());
    cellEdgeCellsPtr_ = new labelListList(nCells());
    
    scalarListList& cellEdgeVols = *cellEdgeVolsPtr_;
    labelListList&  cellEdgeCells = *cellEdgeCellsPtr_;

    // Loop around each primal cell to calculate edgeVols and dual edgeCells
    forAll(cellEdgeVols, cellI)
    {
        scalarList& edgeVols = cellEdgeVols[cellI];
        labelList& edgeCells = cellEdgeCells[cellI];
        
        const cell& c = cells()[cellI];
        // set the number of vertices of the primal cell 2d shape (ie the
        // number of vertical edges)
        const label nVerts = c.size()-2;
        edgeVols.setSize(nVerts);
        edgeCells.setSize(nVerts);
        
        // hash set to store the dual edge cells as we go along
        labelHashSet edgeCellsTmp;
        
        // Also find the radii of the faces above and below cellI
        scalar Rup = -1;
        scalar Rdown = -1;
        
        // Loop around the vertical faces of the cell
        for(label fi = 0; fi < c.size(); fi++)
        {
            label faceI = c[fi];
            label faced = dualFaceMap()[faceI];
            if (faced != -1)
            {
                label own = dualMesh().owner()[faced];
                if (!edgeCellsTmp.found(own)) edgeCellsTmp.insert(own);
                label nei = dualMesh().neighbour()[faced];
                if (!edgeCellsTmp.found(nei)) edgeCellsTmp.insert(nei);
            }
            else if (Rup < 0) Rup = mag(Cf()[faceI]);
            else Rdown = mag(Cf()[faceI]);
        }
        // Radii difference for calculating volume
        const scalar Rdiff = Rup > Rdown ? (pow(Rup,3)-pow(Rdown,3))/3.:
                                           (pow(Rdown,3)-pow(Rup,3))/3.;
        
        // Check that edgeCellsTmp is the correct size
        if (edgeCellsTmp.size() != edgeCells.size())
        {
            FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                << "edgeCellsTmp.size() = " << edgeCellsTmp.size()
                << " edgeCells.size() = " << edgeCells.size()
                << " should be the same size" << abort(FatalError);
        }
        edgeCells = edgeCellsTmp.toc();
        
        // Now loop through the edge cells to set the edgeVols
        for(label iv = 0; iv < edgeCells.size(); iv++)
        {
            label cellv = edgeCells[iv];
            // Find the two faces that cellI and cellv have in common
            // (faceI and faceJ) (correspoding to faceD and faceE on dual)
            label faceI = -1;
            label faceJ = -1;
            const cell& cd = dualMesh().cells()[cellv];
            for(label fd = 0; fd < cd.size(); fd++)
            {
                label faced = cd[fd];
                for(label fi = 0; fi < c.size(); fi++)
                {
                    if (faced == dualFaceMap()[c[fi]])
                    {
                        if (faceI == -1)
                        {
                            faceI = c[fi];
                        }
                        else if(faceJ == -1)
                        {
                            faceJ = c[fi];
                        }
                        else
                        {
                            FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                   << "primal and dual cells should not share more than 2 faces"
                            << abort(FatalError);
                        }
                    }
                }
            }
            // Check that faceI and faceJ are both set
            if (faceI == -1 || faceJ == -1)
            {
                FatalErrorIn("fvMeshWithDual::makeCellEdgeInfo")
                    << "2 faces in common not found for primal cell " << cellI
                    << " and dual cell " << cellv << abort(FatalError);
            }
            
            label faceD = dualFaceMap()[faceI];
            label faceE = dualFaceMap()[faceJ];

            // Find one of the points that faceI and faceJ have in common
            label ip = -1;
            const face& facePointsI = faces()[faceI];
            const face& facePointsJ = faces()[faceJ];
            for (label ipi = 0; ipi < facePointsI.size() && ip == -1; ipi++)
            {
                label fpi = facePointsI[ipi];
                for(label ipj = 0; ipj < facePointsJ.size() && ip == -1; ipj++)
                {
                    if (facePointsJ[ipj] == fpi) ip = fpi;
                }
            }
            
            // Find one of the points that faceD and faceE have in common
            label id = -1;
            const face& facePointsD = dualMesh().faces()[faceD];
            const face& facePointsE = dualMesh().faces()[faceE];
            for (label ipi = 0; ipi < facePointsE.size() && id == -1; ipi++)
            {
                label fpi = facePointsD[ipi];
                for(label ipj = 0; ipj < facePointsE.size() && id == -1; ipj++)
                {
                    if (facePointsE[ipj] == fpi) id = fpi;
                }
            }
            
            // We now have primal cell cellI, dual cell cellv, primal faces
            // faceI and faceJ matching dual faces faceD and faceE
            // From these calculate the edgeVols
            //const point& Cp = C()[cellI];
            //const point& Cd = dualMesh().C()[cellv];
            const point& Cp = dualMesh().points()[id];
            const point& Cd = points()[ip];
            const point& Cf0 = intersections()[faceI];
            const point& Cf1 = intersections()[faceJ];
            edgeVols[iv] = Rdiff*
            (
                sphTriSolidAngle(Cp,Cd,Cf0)
              + sphTriSolidAngle(Cp,Cd,Cf1)
            );
        }
//        // Correct so that edgeVols sum to cell vols
//        scalar sumVol = 0;
//        for(label iv = 0; iv < edgeVols.size(); iv++)
//        {
//            sumVol += edgeVols[iv];
//        }
//        for(label iv = 0; iv < edgeVols.size(); iv++)
//        {
//            edgeVols[iv] *= V()[cellI]/sumVol;
//        }
//        Info << "%age difference between primal vols and sum edge vols = "
//             << 100*(sumVol - V()[cellI])/V()[cellI] << nl;
    }
    
//    // Scale edgeVols so that they sum to dual cell vols
//    scalarList dualVols(dualMesh().nCells(), scalar(0));
//    for(label cellI = 0; cellI < nCells(); cellI++)
//    {
//        for(label iv = 0; iv < cellEdgeVols[cellI].size(); iv++)
//        {
//            dualVols[cellEdgeCells[cellI][iv]] += cellEdgeVols[cellI][iv];
//        }
//    }
//    for(label cellI = 0; cellI < nCells(); cellI++)
//    {
//        for(label iv = 0; iv < cellEdgeVols[cellI].size(); iv++)
//        {
//            label cellv = cellEdgeCells[cellI][iv];
//            cellEdgeVols[cellI][iv] *= dualMesh().V()[cellv]/dualVols[cellv];
//        }
//    }
//    Info << "%age difference between dual vols and dual vol sums = "
//         << (dualVols-dualMesh().V())/dualMesh().V()*100 << endl;
}

void fvMeshWithDual::makerHat() const
{
    if (rHatPtr_) FatalErrorIn("fvMeshWithDual::makerHat")
                      << "rHat already exists" << abort(FatalError);
    rHatPtr_ = new volVectorField("rHat", unitVector(C()));
}

void fvMeshWithDual::makelonHat() const
{
    if (lonHatPtr_) FatalErrorIn("fvMeshWithDual::makelonHat")
                        << "lonHat already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonHatPtr_ = new volVectorField("lonHat", unitVector(k ^ rHat()));
}

void fvMeshWithDual::makelatHat() const
{
    if (latHatPtr_) FatalErrorIn("fvMeshWithDual::makelatHat")
                        << "latHat already exists" << abort(FatalError);
    latHatPtr_ = new volVectorField("latHat", unitVector(rHat() ^ lonHat()));
}

void fvMeshWithDual::makerHatf() const
{
    if (rHatfPtr_) FatalErrorIn("fvMeshWithDual::makerHatf")
                        << "rHatf already exists" << abort(FatalError);
    //rHatfPtr_ = new surfaceVectorField("rHatf", unitVector(intersections()));
    rHatfPtr_ = new surfaceVectorField("rHatf", unitVector(Cf()));
}

void fvMeshWithDual::makelonHatf() const
{
    if (lonHatfPtr_) FatalErrorIn("fvMeshWithDual::makelonHatf")
                        << "lonHatf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonHatfPtr_= new surfaceVectorField("lonHatf",unitVector(k^rHatf()));
}

void fvMeshWithDual::makelatHatf() const
{
    if (latHatfPtr_) FatalErrorIn("fvMeshWithDual::makelatHatf")
                        << "latHatf already exists" << abort(FatalError);
    latHatfPtr_= new surfaceVectorField("latHatf",unitVector(rHatf()^lonHatf()));
}

void fvMeshWithDual::makelon() const
{
    if (lonPtr_) FatalErrorIn("fvMeshWithDual::makelon")
                        << "lon already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonPtr_ = new volScalarField
    (
        "lon",
        atan2
        (
            C() & unitVector(vector(0.,1.,0.) - k.y()*k),
            C() & unitVector(vector(1.,0.,0.) - k.x()*k)
        )
    );
}

void fvMeshWithDual::makelat() const
{
    if (latPtr_) FatalErrorIn("fvMeshWithDual::makelat")
                        << "lat already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    latPtr_ = new volScalarField("lon", asin(rHat() & k));
}

void fvMeshWithDual::makeHeight() const
{
    if (heightPtr_) FatalErrorIn("fvMeshWithDual::makeHeight")
                        << "height already exists" << abort(FatalError);
    heightPtr_ = new volScalarField("height", mag(C()) - earthRadius_);
}

void fvMeshWithDual::makelonf() const
{
    if (lonfPtr_) FatalErrorIn("fvMeshWithDual::makelonf")
                        << "lonf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    lonfPtr_ = new surfaceScalarField
    (
        "lonf",
        atan2
        (
//            intersections() & unitVector(vector(0.,1.,0.) - k.y()*k),
//            intersections() & unitVector(vector(1.,0.,0.) - k.x()*k)
            Cf() & unitVector(vector(0.,1.,0.) - k.y()*k),
            Cf() & unitVector(vector(1.,0.,0.) - k.x()*k)
        )
    );
}

void fvMeshWithDual::makelatf() const
{
    if (latfPtr_) FatalErrorIn("fvMeshWithDual::makelatf")
                        << "latf already exists" << abort(FatalError);
    vector k = mag(OmegaHat_) < SMALL ? vector(0,0,1) : OmegaHat_;
    latfPtr_ = new surfaceScalarField("lonf", asin(rHatf() & k));
}

void fvMeshWithDual::makeHeightf() const
{
    if (heightfPtr_) FatalErrorIn("fvMeshWithDual::makeHeightf")
                        << "heightf already exists" << abort(FatalError);
//    heightfPtr_ = new surfaceScalarField("heightf", mag(intersections()) - earthRadius_);
    heightfPtr_ = new surfaceScalarField("heightf", mag(Cf()) - earthRadius_);
}

void fvMeshWithDual::makeidir() const
{
    if (idirPtr_) FatalErrorIn("fvMeshWithDual::makeidir")
                    << "idir already exists" << abort(FatalError);
    idirPtr_ = new surfaceVectorField("idir", unitVector(Sf()));
}

void fvMeshWithDual::makejdir() const
{
    if (jdirPtr_) FatalErrorIn("fvMeshWithDual::makejdir")
                    << "jdir already exists" << abort(FatalError);
    jdirPtr_ = new surfaceVectorField("jdir", kdir() ^ idir());
    
//    surfaceVectorField& j = *jdirPtr_;
//    
//    forAll(j, faceI)
//    {
//        scalar magj = mag(j[faceI]);
//        if (magj < 0.1) j[faceI] = idir()[faceI];
//        else j[faceI] /= magj;
//    }

//    forAll(j.boundaryField(), patchI)
//    {
//        forAll(j.boundaryField()[patchI], patchFace)
//        {
//            vector& ji = j.boundaryField()[patchI][patchFace];
//            scalar magj = mag(ji);
//            if (magj < 0.1)
//            {
//                ji = idir().boundaryField()[patchI][patchFace];
//            }
//            else ji /= magj;
//        }
//    }
}

void fvMeshWithDual::makekdir() const
{
    if (kdirPtr_) FatalErrorIn("fvMeshWithDual::makekdir")
                    << "kdir already exists" << abort(FatalError);
    if (isSpherical())
    {
        //kdirPtr_ = new surfaceVectorField("kdir", unitVector(intersections()));
        kdirPtr_ = new surfaceVectorField("kdir", unitVector(Cf()));

        surfaceVectorField& k = *kdirPtr_;
        
        forAll(k, faceI)
        {
            if (faceToPatchEdge()[faceI] == -1) // horizontal face
            k[faceI] = vector(k[faceI].y(), k[faceI].z(), k[faceI].x());
        }

        forAll(k.boundaryField(), patchI)
        {
            forAll(k.boundaryField()[patchI], patchFace)
            {
                vector& ki = k.boundaryField()[patchI][patchFace];
                if
                (
                    faceToPatchEdge()[boundaryMesh()[patchI].start()+patchFace]
                 == -1
                )
                {
                    ki = vector(ki.y(), ki.z(), ki.x());
                }
            }
        }
    }
    else // not spherical
    {
        kdirPtr_ = new surfaceVectorField
        (
            IOobject
            (
                "kdir", time().findInstance(meshDir(), "points"),
                meshSubDir, *this
            ),
            *this,
            dimensionedVector("kdir", dimless, OmegaHat())
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
