/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::dualDecomp

Description
    For dual meshes. Follows primal decomposition - gets point decomposition
    from the primal decomposition and then decomposes dual according to minimum
    point decomposition. Needs
    - all extrusion maps from cell/face/point back to patch (see meshStructure
      application
    - constant/cellDecomposition as written by 'decomposePar -cellDist'

    Run on the dual mesh, specify the primal region and patch.

SourceFiles
    dualDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef dualDecomp_H
#define dualDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dualDecomp Declaration
\*---------------------------------------------------------------------------*/

class dualDecomp
:
    public decompositionMethod
{
    // Private data

        dictionary methodDict_;

        word region_;

        word patch_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const dualDecomp&);
        dualDecomp(const dualDecomp&);


public:

    //- Runtime type information
    TypeName("dual");


    // Constructors

        //- Construct given the decomposition dictionary
        dualDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~dualDecomp()
    {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
