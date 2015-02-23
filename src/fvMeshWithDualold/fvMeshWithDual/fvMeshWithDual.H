/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open So
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

Class
    Foam::fvMeshWithDual

Description
    Re-reads the face and cell centres, areas and volumes or calculates them
    using spherical geometry and also reads addressing from each
    cell and face onto a patch

SourceFiles
    fvMeshWithDual.C
    fvMeshWithDualGeometry.C
    fvMeshWithDualnonOrthog.C
    fvMeshWithDualTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshWithDual_H
#define fvMeshWithDual_H

#include "primitiveMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvMeshWithDual Declaration
\*---------------------------------------------------------------------------*/

class fvMeshWithDual
:
    public fvMesh
{
    // Private data
    
        // The earthProperties dictionary which holds Omega and earthRadius
        IOdictionary earthProperties_;
    
        // Is this a spherical mesh or not?
        const Switch isSpherical_;
    
        //- name of the patch for the surface of the sphere
        const word patchName_;
        
        //- patch for the surface of the sphere
        const polyPatch& bottomPatch_;
        
        // Rotation of the sphere
        const dimensionedVector Omega_;
        
        // Magnitude of gravity 
        const dimensionedScalar magg_;
    
        // The axis of rotation of the sphere
        const vector OmegaHat_;
        
        // The reference height
        const dimensionedScalar earthRadius_;
        
        //- Number of levels (or layers) of mesh
        const label nLevels_;

        //- Map from each vertical face onto a patch edge
        mutable IOList<label>* faceToPatchEdgePtr_;

        //- Map from each horizontal face onto a patch face
        mutable IOList<label>* faceToPatchFacePtr_;

        //- Map from each point onto a patch point
        mutable IOList<label>* pointToPatchPointPtr_;
        
        //- Map from each cell onto a patch face
        mutable IOList<label>* cellToPatchFacePtr_;

        //- Reference to the dual mesh (ie the dual of this mesh so it could
        // be the primal of the simulation mesh)
        const fvMeshWithDual* dualMesh_;

        //- Map from each (vertical) face to its dual
        mutable IOList<label>* dualFaceMapPtr_;
        
        //- Orientation of face map from primal to dual
        mutable List<bool>* signFaceMapPtr_;
        
        //- volumes associated with each face corresponding to the owner
        //  and neighbour cell
        mutable surfaceScalarField* faceVolOwnPtr_;
        mutable surfaceScalarField* faceVolNeiPtr_;
        mutable surfaceScalarField* faceVolPtr_;
        
        //- height (depth) of each vertical face (in the vertical)
        mutable surfaceScalarField* depthfPtr_;
        
        //- area fractions associated with the bottom patch
        mutable List<List<FixedList<scalar,2> > >* fpeAreaFracPtr_;
        mutable scalarListList* facePointAreaFracPtr_;
        
        //- list of points of each face of the patch in same order as above
        mutable labelListList* facePointsPtr_;
        
        //- areas associated with each edge of the patch on each side
        mutable scalarList* edgeAreaOwnPtr_;
        mutable scalarList* edgeAreaNeiPtr_;
        
        //- The diagonal of the H matrix for converting fluxes between primal
        //  and dual on non-orthonogal meshes
        mutable surfaceScalarField* HdiagPtr_;
        
        // The intersection points between primal and dual (near face centres)
        mutable surfaceVectorField* intersectionsPtr_;

        // The volumes associated with each vertical edge of each cell
        mutable scalarListList* cellEdgeVolsPtr_;
        
        // The dual cells associated with each vertical edge of each primal cell
        mutable labelListList* cellEdgeCellsPtr_;

        // Unit vectors in the lon, lat and height directions
        mutable volVectorField* rHatPtr_;
        mutable volVectorField* lonHatPtr_;
        mutable volVectorField* latHatPtr_;
        mutable surfaceVectorField* rHatfPtr_;
        mutable surfaceVectorField* lonHatfPtr_;
        mutable surfaceVectorField* latHatfPtr_;

        // The latitude, longitude and height above reference for each
        // cell and face of the mesh
        mutable volScalarField* lonPtr_;
        mutable volScalarField* latPtr_;
        mutable volScalarField* heightPtr_;
        mutable surfaceScalarField* lonfPtr_;
        mutable surfaceScalarField* latfPtr_;
        mutable surfaceScalarField* heightfPtr_;
        
        // i, j, k dir for each vertical face and ddir (cell centre to centre)
        mutable surfaceVectorField* idirPtr_;
        mutable surfaceVectorField* jdirPtr_;
        mutable surfaceVectorField* kdirPtr_;
        mutable surfaceVectorField* ddirPtr_;
        
    // Private Member Functions

        //- Disallow construct as copy
        fvMeshWithDual(const fvMeshWithDual&);

        //- Disallow assignment
        void operator=(const fvMeshWithDual&);

        //- Read in new face centres and areas if present or calculate
        //- based on spherical geometry
        void newFaceCentresAndAreas(const bool readOverrides);

        //- Read in new cell centres and volumes if present or calculate
        //- based on spherical geometry
        void newCellCentresAndVols(const bool readOverrides);
        
        
        // Read demand driven data
            void readFaceToPatchEdge() const;
            void readFaceToPatchFace() const;
            void readPointToPatchPoint() const;
            void readCellToPatchFace() const;
            void readDualFaceMap() const;
        
        // Make geometric data
            void makesignFaceMap() const;
            void makeSphFaceVols() const;
            void makeFaceDepth() const;
            void makeSphAreaFracs() const;
            void makeHdiag() const;
            void makeIntersections() const;
            void makeCellEdgeInfo() const;
            void makerHat() const;
            void makelonHat() const;
            void makelatHat() const;
            void makerHatf() const;
            void makelonHatf() const;
            void makelatHatf() const;
            void makelon() const;
            void makelat() const;
            void makeHeight() const;
            void makelonf() const;
            void makelatf() const;
            void makeHeightf() const;
            void makeidir() const;
            void makejdir() const;
            void makekdir() const;
            void makeddir() const;

        // Helper functions
        const polyPatch& findBottomPatch(const word patchName) const;

public:

    // Static data members

    // Declare name of the class and its debug switch
    ClassName("fvMeshWithDual");

    // Constructors

        //- Construct from IOobject
        explicit fvMeshWithDual
        (
            const IOobject& io,
            const bool readOverrides = true
        );

        //- Construct from IOobject and reference to the dual
        explicit fvMeshWithDual
        (
            const IOobject& io,
            const fvMeshWithDual& dualMesh,
            const bool readOverrides = true
        );


    // Member Functions

        // Access

            // Is this a spherical mesh or not?
            Switch isSpherical() const { return isSpherical_;}

            // The earthProperties dictionary
            const IOdictionary& earthProperties() const
                { return earthProperties_;}
        
            //- name of the patch for the surface of the sphere
            const word& patchName() const { return patchName_; }

            //- patch for the surface of the sphere
            const polyPatch& bottomPatch() const {return bottomPatch_; }

            // Rotation of the sphere
            const dimensionedVector& Omega() const { return Omega_;}
            
            // Magnitude of gravity 
            const dimensionedScalar& magg() const { return magg_;}

            // The axis of rotation of the sphere
            const vector& OmegaHat() const { return OmegaHat_;}
        
            // The reference height
            const dimensionedScalar& earthRadius() const { return earthRadius_;}

            // Number of levels (or layers) of mesh
            label nLevels() const { return nLevels_; }

            //- Map from each vertical face onto a patch edge
            const IOList<label>& faceToPatchEdge() const
            {
                if(!faceToPatchEdgePtr_) readFaceToPatchEdge();
                return *faceToPatchEdgePtr_;
            }

            //- Map from each horizontal face onto a patch face
            const IOList<label>& faceToPatchFace() const
            {
                if (!faceToPatchFacePtr_) readFaceToPatchFace();
                return *faceToPatchFacePtr_;
            }

            //- Map from each point onto a patch point
            const IOList<label>& pointToPatchPoint() const
            {
                if (!pointToPatchPointPtr_) readPointToPatchPoint();
                return *pointToPatchPointPtr_;
            }

            //- Map from each cell onto a patch face
            IOList<label> cellToPatchFace() const
            {
                if (!cellToPatchFacePtr_) readCellToPatchFace();
                return *cellToPatchFacePtr_;
            }

            //- The dual of this mesh
            const fvMeshWithDual& dualMesh() const
            {
                if (!dualMesh_)
                {
                    FatalErrorIn("fvMeshWithDual::dualMesh")
                        << "dual mesh not set" << abort(FatalError);
                }
                return *dualMesh_;
            }

            //- Set the dual of this mesh (modifies face areas of dual for
            //  consistency with primal
            void setDual(const fvMeshWithDual& dm);

            //- Map from each (vertical) face to its dual (-1 if horizontal face)
            const IOList<label>& dualFaceMap() const
            {
                if (!dualFaceMapPtr_) readDualFaceMap();
                return *dualFaceMapPtr_;
            }

            //- Orientation of face map from primal to dual
            const List<bool>& signFaceMap() const
              {if(!signFaceMapPtr_) makesignFaceMap(); return *signFaceMapPtr_;}

            //- The sign of the orientation from primal to dual for faceI
            int signMap(const label faceI) const
                {return 2*signFaceMap()[faceI]-1; }

            //- volumes associated with each face corresponding to the owner
            const surfaceScalarField& faceVolOwn() const
                {if(!faceVolOwnPtr_) makeSphFaceVols(); return *faceVolOwnPtr_;}

            //- volumes associated with each face corresponding to the neighbour
            const surfaceScalarField& faceVolNei() const
                {if(!faceVolNeiPtr_) makeSphFaceVols(); return *faceVolNeiPtr_;}

            //- total volumes associated with each face
            const surfaceScalarField& faceVol() const
                {if (!faceVolPtr_) makeSphFaceVols(); return *faceVolPtr_;}
                
            //- height (depth) of each vertical face (in the vertical)
            const surfaceScalarField& depthf() const
                {if (!depthfPtr_) makeFaceDepth(); return *depthfPtr_;}
            
            //- area fraction of each face of the patch associated with each
            //  patch for each edge
            const List<List<FixedList<scalar,2> > >& fpeAreaFrac() const
               {if(!fpeAreaFracPtr_)makeSphAreaFracs();return *fpeAreaFracPtr_;}

            //- area fraction of each face of patch associated with each point
            const scalarListList& facePointAreaFrac() const
   {if(!facePointAreaFracPtr_)makeSphAreaFracs();return *facePointAreaFracPtr_;}

            //- points of faces of patch in the same order as facePointAreaFrac
            const labelListList& facePoints() const
                {if(!facePointsPtr_)makeSphAreaFracs(); return *facePointsPtr_;}

            //- areas associated with each edge of the patch on owner side
            const scalarList& edgeAreaOwn() const
               {if(!edgeAreaOwnPtr_)makeSphAreaFracs();return *edgeAreaOwnPtr_;}

            //- areas associated with each edge of the patch on neighbour side
            const scalarList& edgeAreaNei() const
               {if(!edgeAreaNeiPtr_)makeSphAreaFracs();return *edgeAreaNeiPtr_;}

            //- The diagonal of the H matrix for converting fluxes between
            //  primal and dual on non-orthonogal meshes
            const surfaceScalarField& Hdiag() const
                { if (!HdiagPtr_) makeHdiag(); return *HdiagPtr_;}

            //- The intersection points between primal and dual (near faceCentres)
            const surfaceVectorField& intersections() const
            {
                if(!intersectionsPtr_) makeIntersections();
                return *intersectionsPtr_;
            }

            //- The volumes associated with each vertical edge of each cell
            const scalarListList& cellEdgeVols() const
            {
                if (!cellEdgeVolsPtr_) makeCellEdgeInfo();
                return *cellEdgeVolsPtr_;
            }
            
            // Dual cells associated with each vertical edge of each primal cell
            const labelListList& cellEdgeCells() const
            {
                if (!cellEdgeCellsPtr_) makeCellEdgeInfo();
                return *cellEdgeCellsPtr_;
            }

            // The latitude, longitude and height above reference for each
            // cell and face of the mesh
            const volScalarField& lon() const
                {if (!lonPtr_) makelon(); return *lonPtr_;}
            const volScalarField& lat() const
                { if (!latPtr_) makelat(); return *latPtr_;}
            const volScalarField& height() const
                { if (!heightPtr_) makeHeight(); return *heightPtr_;}
            const surfaceScalarField& lonf() const
                { if (!lonfPtr_) makelonf(); return *lonfPtr_;}
            const surfaceScalarField& latf() const
                { if (!latfPtr_) makelatf(); return *latfPtr_;}
            const surfaceScalarField& heightf() const
                { if (!heightfPtr_) makeHeightf(); return *heightfPtr_;}
            
            // Unit vectors in the lon, lat and height directions
            const volVectorField& lonHat() const
                { if (!lonHatPtr_) makelonHat(); return *lonHatPtr_; }
            const volVectorField& latHat() const
                { if (!latHatPtr_) makelatHat(); return *latHatPtr_; }
            const volVectorField& rHat() const
                { if (!rHatPtr_) makerHat(); return *rHatPtr_; }
            const surfaceVectorField& lonHatf() const
                { if (!lonHatfPtr_) makelonHatf(); return *lonHatfPtr_; }
            const surfaceVectorField& latHatf() const
                { if (!latHatfPtr_) makelatHatf(); return *latHatfPtr_; }
            const surfaceVectorField& rHatf() const
                { if (!rHatfPtr_) makerHatf(); return *rHatfPtr_; }
            const surfaceVectorField& idir() const
                { if (!idirPtr_) makeidir(); return *idirPtr_; }
            const surfaceVectorField& jdir() const
                { if (!jdirPtr_) makejdir(); return *jdirPtr_; }
            const surfaceVectorField& kdir() const
                { if (!kdirPtr_) makekdir(); return *kdirPtr_; }
            const surfaceVectorField& ddir() const
                { if (!ddirPtr_) makeddir(); return *ddirPtr_; }

            //- Return face deltas as surfaceSclarField
            tmp<surfaceScalarField> magDelta() const;

        // Map fields from one mesh to another
        
            //- Map surfaceField from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;

            //- Map signed surfaceField from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > signedDualMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > signedDualMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;

            //- Map surfaceField of fluxes from this mesh to its dual
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualFluxMap
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > dualFluxMap
            (
                const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMeshWithDualTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
