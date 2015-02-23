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
    Foam::TRiSK3dData

Description
    All the data necessary for TRiSK3d operators in 3D

SourceFiles
    TRiSK3dData.C

\*---------------------------------------------------------------------------*/

#ifndef TRiSK3dData_H
#define TRiSK3dData_H

#include "MeshObject.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "pointFields.H"
#include "sphericalGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class TRiSK3dData Declaration
\*---------------------------------------------------------------------------*/

class TRiSK3dData
:
    public MeshObject<fvMesh, MoveableMeshObject, TRiSK3dData>
{
    // Private data

        //- The unit patch on which geometry is calculated;
        PrimitivePatch<face, List, pointField> unitPatch_;

        //- Number of mesh levels
        const label nLevels_;
        
        //- For each cell of the mesh, the k dir
        const volVectorField rHat_;

        //- For each cell of the mesh, the k direction
        const surfaceVectorField rHatf_;
                
        //- For each edge of face, the normal direction
        const surfaceVectorField idir_;
        
        //- For each face, the tangential direction
        const surfaceVectorField jdir_;

        //- The lat and lon directions for each cell and each face
        const volVectorField lonHat_;
        const volVectorField latHat_;
        const surfaceVectorField lonHatf_;
        const surfaceVectorField latHatf_;

        //- The areas on each side of each face and their sum
        surfaceScalarField faceVolOwn_;
        surfaceScalarField faceVolNei_;
        surfaceScalarField faceVol_;
                

    // Private Member Functions

        //- Calculate areas and area fractions, volumes and depths
        void calcAreasVols();
        
        //- Calculate the delta coeffs
        void calcDeltaCoeffs();

        //- Check the delta coeffs
        void checkDeltaCoeffs();

public:

    // Static data members

        //- name of the patch from which the unit patch is constructed
        static word patchName;
        
    // Declare name of the class and its debug switch
    TypeName("TRiSK3dData");

    // Constructors
    
        //- Construct from mesh
        explicit TRiSK3dData(const fvMesh& mesh);

    //- Destructor
    virtual ~TRiSK3dData()
    {}

    // Member Functions

        // Access

            //- The unit patch on which geometry is calculated;
            const PrimitivePatch<face, List, pointField>& unitPatch() const
            {
                return unitPatch_;
            }

            //- Number of mesh levels
            label nLevels() const { return nLevels_; }
            
            //- For each cell of the mesh, the k dir
            const volVectorField& rHat() const { return rHat_; }

            //- For each face of the mesh, the k direction
            const surfaceVectorField& rHatf() const { return rHatf_; }
            
            //- For each face, the normal direction
            const surfaceVectorField& idir() const { return idir_; }
            
            //- For each face, the tangential direction
            const surfaceVectorField& jdir() const { return jdir_; }
            
            //- the lat and lon directions for each cell and face
            const volVectorField& lonHat() const { return lonHat_; }
            const volVectorField& latHat() const { return latHat_; }
            const surfaceVectorField& lonHatf() const { return lonHatf_; }
            const surfaceVectorField& latHatf() const { return latHatf_; }
                        
            //- The areas on each side of each face and their sum
            const surfaceScalarField& faceVolOwn() const { return faceVolOwn_;}
            
            const surfaceScalarField& faceVolNei() const { return faceVolNei_;}

            const surfaceScalarField& faceVol() const { return faceVol_;}

        //- Delete the data when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
