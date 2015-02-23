/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::meshToPointField

Description
    Abstract base type for holding the stencils of cells and
    interpolation weights for interpolating from the cell centres of
    a mesh onto a field of points

SourceFiles
    meshToPointField.C
    meshToPoinI.H

\*---------------------------------------------------------------------------*/

#ifndef meshToPointField_H
#define meshToPointField_H

#include "meshToPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToPointField Declaration
\*---------------------------------------------------------------------------*/

template<class approxType>
class meshToPointField
{
    // Private data

        const fvMesh& mesh_;
        
        labelListList stencils_;
        
        scalarListList weights_;
        
    public:
        enum stencilType{CELLCELLS, CELLPOINTCELLS, UPWINDCPC};
    private:
        const stencilType stencilType_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshToPointField(const meshToPointField<approxType>&);

        //- Disallow default bitwise assignment
        void operator=(const meshToPointField<approxType>&);


public:

    // Constructors

        //- Construct from components, specifying type of stencil
        meshToPointField
        (
            const pointField& pts,
            const fvMesh& mesh,
            const stencilType sType = CELLCELLS
        );

        //- Construct from mesh (need to set point before using)
        meshToPointField
        (
            const fvMesh& mesh,
            const stencilType sType = CELLCELLS
        );


    //- Destructor
    ~meshToPointField();


    // Access functions
    
        //- interpolation stencils
        const labelListList& stencils() const { return stencils_; }

        //- interpolation weights
        const scalarListList& weights() const { return weights_; }

    // Member Functions

        //- Interpolate volField onto these points
        template<class Type>
        inline Field<Type> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
        //- Calculate the stencil and weights for CELLCELLS stencils
        void setPointsCellCells(const pointField& pts);
        
        //- Calculate the stencil and weights for CELLPOINTCELLS stencils
        void setPointsCellPointCells(const pointField& pts);

        //- Calculate the stencil and weights for UPWINDCPC stencils
        void setPointsUpwindCPC
        (
            const pointField& pts
        );

        //- Calculate the stencil and weights based on new points
        void setPoints(const pointField& pts)
        {
            if (stencilType_ == CELLCELLS) setPointsCellCells(pts);
            else if (stencilType_ == CELLPOINTCELLS)
                setPointsCellPointCells(pts);
            else if (stencilType_ == UPWINDCPC) setPointsUpwindCPC(pts);
            else
            {
                FatalErrorIn("meshToPointField::setPoints") << stencilType_
                    << " is not a stencil type" << endl;
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToPointFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshToPointField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
