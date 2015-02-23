/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::TRiSKData

Description
    All the data necessary for TRiSK operators to be used with fvMeshWithDual

SourceFiles
    TRiSKData.C
    TRiSKcircToFluxDiagData.C

\*---------------------------------------------------------------------------*/

#ifndef TRiSKData_H
#define TRiSKData_H

#include "meshWithDual.H"
#include "MeshObject.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class TRiSKData Declaration
\*---------------------------------------------------------------------------*/

class TRiSKData
:
    public MeshObject<fvMeshWithDual, MoveableMeshObject, TRiSKData>
{
    // Private data

        const fvMeshWithDual& smesh_;
        
    public:
        enum nonOrthogHtype{DIAGONAL, DUBOS, DUBOSW, ASYMMETRIC,NONE};
    private:        
        nonOrthogHtype Htype_;

//        //- The solid angle subtended my each column (face areas of bottomPatch)
//        scalarField solidAngle_;
//        
//        //- The ratio of each patch face area associated with each point
//        scalarListList facePointAreaFracs_;
//        
//        //- Area fraction of each dual face associated with each primal face
//        scalarListList dualAreaFracs_;

//        //- The points of each patch face (same order as facePointAreaFracs)
//        labelListList facePoints_;

        //- The stencil of faces in the owner cell for each face of the mesh
        labelListList ownerStencil_;

        //- The stencil of faces in the neighbour cell for each face of the mesh
        labelListList neighbourStencil_;

        //- The stencil of horizontal faces to map velocity to vertical faces
        List<FixedList<label, 4> > vertToHorizStencil_;
        labelList vertToHorizStencilSize_;

        //- For each face on the mesh store the values which multiply the
        //  normal fluxes of the owner cell faces to obtain the perpendicular
        scalarListList own_weights_;
        
        //- For each face on the mesh store the values which multiply the
        //  normal fluxes of the neighbour cell faces to obtain the perpendicular
        scalarListList nei_weights_;
        
        //- The diagonal entries in the H operator for converting from dual
        //- circulations to primal fluxes
        surfaceScalarField circToFluxDiag_;
        
        //- The off-diagonal entries in the H operator for converting from dual
        //- circulations to primal fluxes
        List<FixedList<scalar, 6> > circToFluxOffDiagWeights_;
        
        //- The stencil of each face for circToFluxOffDiag
        List<FixedList<label, 6> > circToFluxStencil_;
        labelList circToFluxStencilSize_;
        
        //- The off-diagonal entries in the H operator for converting from
        //- primal velocity in ddir to idir
        labelListList ddirToFluxOffDiagWeights_;
        
        //- The stencil of each face for ddirToFluxOffDiag
        labelListList ddirToFluxStencil_;
        
    // Private Member Functions

//        //- Calculate areas and area fractions, volumes
//        void calcAreasVols();
//        
//        //- Calculate the delta coeffs for primal and dual meshes
//        void calcDeltaCoeffs();

//        //- Check the delta coeffs for primal and dual meshes
//        void checkDeltaCoeffs();

        //- Compare primal and dual mesh points and cell centres
        void checkPrimalDual();

        //- Find the face stencil and the weights for face ie
        void calcStencilWeights
        (
            labelList& stencil,
            scalarList& weight,
            const label faceI,
            const label cellI
        ) const;
        
        //- Calculate the length, ell, and distance, d of edge ie
        void calcEllDist(scalar& ell, scalar& dist, const label ie) const;
        
        //- Calculate the H operator stencil and weights (the operator for
        //- converting between dual circulations to primal fluxes)
        void calcCircToFluxWeights();
        
        //- Different methods for calculating the H operator
        // H operator from Thomas Dubos
        void calcDubosH
        (
            FixedList<scalar,6>& H,
            FixedList<label,6>& stencil,
            scalar& Hdiag,
            const label faceI
        );
        
        // H operator from Thomas Dubos modified by Weller
        void calcDubosWH
        (
            FixedList<scalar,6>& H,
            FixedList<label,6>& stencil,
            scalar& Hdiag,
            const label faceI
        );
        
        // Reconstruction exact at the centre
        void reconstructH
        (
            FixedList<scalar,6>& H,
            FixedList<label,6>& stencil,
            scalar& Hdiag,
            const label faceI,
            const bool compact
        );
        
        // Reconstruction exact at the centre from primal only (3d)
        void reconstructH
        (
            labelList& Hp,
            labelList& stencilp,
            scalar& Hdiag,
            const label faceI
        );
        
        // Make H symmetric
        void symmetriseH();

        // The owner cell (face) of patch edge ie
        label edgeOwner(const label ie) const
            { return smesh().bottomPatch().edgeFaces()[ie][0]; }
        label edgeNeib(const label ie)const
            { return smesh().bottomPatch().edgeFaces()[ie].last(); }


public:

    // Declare name of the class and its debug switch
    TypeName("TRiSKData");


    // Constructors
    
        //- Construct from mesh
        explicit TRiSKData(const fvMesh& mesh);

    //- Destructor
    virtual ~TRiSKData()
    {}

    // Member Functions

        // Access
        
            //- The fvMeshWithDual for this field
            const fvMeshWithDual& smesh() const { return smesh_; }
            
            //- The type of non-orthogonal corrector, H
            const nonOrthogHtype& Htype() const { return Htype_; }

//            //- Solid angle subtended my each column (face areas of bottomPatch)
//            const scalarField& solidAngle() const { return solidAngle_; }
//            
//            //- The ratio of each patch face area associated with each point
//            const scalarListList& facePointAreaFracs() const
//            {
//                return facePointAreaFracs_;
//            }
//            
//            //- Area fraction of each dual face associated with each primal face
//            const scalarListList& dualAreaFracs() const
//            {
//                return dualAreaFracs_;
//            }

//            //- The points of each patch face (same order as facePointAreaFracs)
//            const labelListList& facePoints() const { return facePoints_; }

            //- The stencil of faces in the owner cell for each face of the mesh
            const labelListList& ownerStencil() const {return ownerStencil_;}

            //- The faces in the neighbour cell for each face of the mesh
            const labelListList& neighbourStencil() const
                {return neighbourStencil_;}

            //- The stencil of horizontal faces to map velocity to vertical faces
            const List<FixedList<label, 4> >& vertToHorizStencil() const
                { return vertToHorizStencil_;}
                
            const labelList& vertToHorizStencilSize() const
                { return vertToHorizStencilSize_; }

            //- For each face on the mesh store the values which multiply the
            //  normal fluxes of the owner cell faces to obtain the perpendicular
            const scalarListList& own_weights() const { return own_weights_;}
            
            //- For each face on the mesh store the values which multiply the
            //  fluxes of the neighbour cell faces to obtain the perpendicular
            const scalarListList& nei_weights() const { return nei_weights_;}

            //- The diagonal entries in the H operator for converting from dual
            //- circulations to primal fluxes
            const surfaceScalarField& circToFluxDiag() const
                { return circToFluxDiag_;}

            //- The off-diagonal entries in the H operator for converting from dual
            //- circulations to primal fluxes
            const List<FixedList<scalar, 6> >& circToFluxOffDiagWeights() const
                { return circToFluxOffDiagWeights_;}

            //- The stencil of each face for circToFluxOffDiag
            const List<FixedList<label, 6> >& circToFluxStencil() const
                { return circToFluxStencil_;}
            
            //- The stencil size for circToFluxOffDiag for each face
            const labelList& circToFluxStencilSize() const
                { return circToFluxStencilSize_;}

            //- The off-diagonal entries in the H operator for converting
            //- from primal velocity in ddir to idir
            const labelListList& ddirToFluxOffDiagWeights() const
                { return ddirToFluxOffDiagWeights_;}
            
            //- The stencil of each face for ddirToFluxOffDiag
            const labelListList& ddirToFluxStencil() const
                { return ddirToFluxStencil_;}

        //- Delete the data when the mesh moves
        virtual bool movePoints();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
