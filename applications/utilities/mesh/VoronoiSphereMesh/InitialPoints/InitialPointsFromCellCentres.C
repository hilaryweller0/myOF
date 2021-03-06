/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "InitialPointsFromCellCentres.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "findInterpWeights.H"
//#include "meshToPoint.H"
//#include "polyFit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(InitialPointsFromCellCentres, 0);
addToRunTimeSelectionTable(InitialPointsFrom, InitialPointsFromCellCentres, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

InitialPointsFromCellCentres::InitialPointsFromCellCentres
(
    const IOdictionary& dict
)
:
    InitialPointsFrom(dict),
    case_(dict.lookup("sourceCase")),
    runTime_(Time::controlDictName, fileName(case_.path()), fileName(case_.name())),
    mesh_
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime_.timeName(),
            runTime_,
            IOobject::MUST_READ
        )
    ),
    uniform_(dict.found("uniformResolution")),
    uniformRes_(dict.lookupOrDefault<scalar>("uniformResolution", scalar(1))),
    reqRes_
    (
        IOobject
        (
            "requiredResolution", runTime_.constant(), mesh_,
            IOobject::MUST_READ
        ),
        mesh_
//        volScalarField
//        (
//            IOobject("requiredResolution", runTime_.timeName(), mesh_),
//            mesh_,
//            dimensionedScalar("rr", dimLength, uniformRes_)
//        )
    )
{
    if (uniform_)
    {
        Info << "Constructing InitialPointsFromCellCentres with uniform required resolution " << uniformRes_ << endl;
    }
    else
    {
        Info << "Constructing InitialPointsFromCellCentres with non-uniform required resolution read from " << reqRes_.path()/reqRes_.name() << endl;
    }
//    Info << "reqRes = " << reqRes_ << nl;
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

scalar InitialPointsFromCellCentres::interpolate
(
    const volScalarField& f,
    const point& p
) const
{
//    meshToPoint<polyFit<1> > mtp(p, f.mesh());
//    scalar fp = mtp.interpolate(f);

/*// Tanh function
    point centreHat(0,0.866025404, .5);
    scalar alpha = 0.157079633;
    scalar beta = 0.523598776;
    scalar gamma = 0.000244141;
    scalar r = acos(p & centreHat);
    return 1/pow(0.5/(1+gamma)*(tanh((beta-r)/alpha)+1) + gamma, 0.25);
*/
 // Quadratic fit (clipped)
    label celli = mesh_.findNearestCell(p);
    const labelList& cells = mesh_.cellCells()[celli];
    const pointField& C = mesh_.cellCentres();

    // linear interpolation
    scalarList weights;
    vectorField points(cells.size()+1);
    points[0] = C[celli];
    forAll(cells, cellj)
    {
        points[cellj+1] = C[cells[cellj]];
    }
    calculatePolyWeights(weights, points, p, 2, 2, 1, 1000.);

    scalar fp = weights[0]*f[celli];
    forAll(cells, cellj)
    {
        fp += weights[cellj+1]*f[cells[cellj]];
    }

    // clip unboundedness
    scalar fmin = f[celli];
    scalar fmax = f[celli];
    forAll(cells, cellj)
    {
        fmin = min(fmin, f[cells[cellj]]);
        fmax = max(fmax, f[cells[cellj]]);
    }
    fp = min(max(fp, fmin), fmax);

    return fp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
