/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "VoronoiSphereMeshing.H"
//#include "fvMeshWithCentres.H"
//#include "fvCFD.H"
//#include "meshTools.H"
#include "extrudedMesh.H"
#include "linearRadial.H"
//#include "polyTopoChange.H"
//#include "polyTopoChanger.H"
//#include "edgeCollapser.H"
//#include "IFstream.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField Foam::VoronoiSphereMeshing::dualPoints() const
{
    pointField dPoints(number_of_finite_cells());

    label dualVerti = 0;
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        cit++
    )
    {
        if (sphericalCell(cit))
        {
            Foam::point v = circumcenter(cit);
            v *= radius/mag(v);
            if (cit->cellIndex() <= 0)
            {
                cit->cellIndex() = dualVerti;
                dualVerti++;
            }
            dPoints[cit->cellIndex()] = v;
        }
    }
    //dPoints.resize(dualVerti);
    return dPoints;
}


const Foam::pointField Foam::VoronoiSphereMeshing::triPoints() const
{
    pointField tPoints(number_of_vertices()-1);
    Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    //label i = 0;
    for(vit++; vit != finite_vertices_end(); vit++)
    {
        tPoints[vit->index()-1] = topoint(vit->point());
    }
    return tPoints;
}


const Foam::faceList Foam::VoronoiSphereMeshing::dualFaces() const
{
    // First set the cell numbers if needed
    if (finite_cells_begin()->cellIndex() <= 0)
    {
        label dualVerti = 0;
        for
        (
            Triangulation::Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            cit++
        )
        {
            if (sphericalCell(cit)) cit->cellIndex() = dualVerti++;
        }
    }

    // then define the dual faces
    faceList dFaces(number_of_vertices()-1);
    //label iv = 0;

    Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    vit++;
    do
    {
        // Circulate around the cells incident to this edge (non infinite)
        Cell_circulator ccirc = incident_cells(edgeBetween(vit,vCentral()));
        Cell_circulator cStart = ccirc;

        label nCells = 0;
        label iv = vit->index()-1;

        do
        {
            if (sphericalCell(ccirc))
            {
                label celli = ccirc->cellIndex();
                if (celli >= 0) nCells++;
            }
        } while (--ccirc != cStart);

        dFaces[iv].setSize(nCells);
        label ic = 0;
        do
        {
            if (sphericalCell(ccirc))
            {
                label celli = ccirc->cellIndex();
                if (celli >= 0)
                {
                    dFaces[iv][ic++] = celli;
                }
            }
        } while (--ccirc != cStart);

        //iv++;
    } while (++vit != finite_vertices_end());

    return dFaces;
}

const Foam::faceList Foam::VoronoiSphereMeshing::triFaces() const
{
    faceList tFaces(number_of_finite_cells(), face(3));
    //label ic = 0;
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        FixedList<Foam::point,3> pts;
        FixedList<label,3> fac;
        label ip = 0;

        for(label i = 0; i < 4; i++)
        {
            label cvi = cit->vertex(i)->index();
            if (cvi != 0)
            {
                if (ip >= 3)
                {
                    FatalErrorIn("VoronoiSphereMeshing::triFaces")
                    << "cell does not have exatly 3 vertices on the surface of "
                        << "the sphere.\nVertices are at locations:\n"
                        << topoint(cit->vertex(0)->point()) << " with mag "
                        << mag(topoint(cit->vertex(0)->point())) << nl
                        << topoint(cit->vertex(1)->point()) << " with mag "
                        << mag(topoint(cit->vertex(1)->point())) << nl
                        << topoint(cit->vertex(2)->point()) << " with mag "
                        << mag(topoint(cit->vertex(2)->point())) << nl
                        << topoint(cit->vertex(3)->point()) << " with mag "
                        << mag(topoint(cit->vertex(3)->point())) << nl
                        << exit(FatalError);
                }

                fac[ip] = cvi-1;
                pts[ip] = topoint(cit->vertex(i)->point());
                ip++;
            }
        }
        label ic = cit->cellIndex();
        if ((pts[0] & (pts[1] ^ pts[2])) < 0)
        {
            tFaces[ic][0] = fac[0];
            tFaces[ic][1] = fac[2];
            tFaces[ic++][2] = fac[1];
        }
        else
        {
            tFaces[ic][0] = fac[0];
            tFaces[ic][1] = fac[1];
            tFaces[ic++][2] = fac[2];
        }
    }
    return tFaces;
}


void Foam::VoronoiSphereMeshing::write
(
    const Time& runTime,
    const word& regionName,
    const bool writeDelaunay
) const
{
    Info << "Writing the primal and dual meshes with centres" << endl;
    // Create the primal mesh patch to extrude
    PrimitivePatch<face, List, pointField> newPatch = !writeDelaunay ?
    PrimitivePatch<face, List, pointField>(dualFaces(), dualPoints()) :
    PrimitivePatch<face, List, pointField>(triFaces(), triPoints());
    
//    // If we actually want the Delaunay mesh, make it the dual
//    if (writeDelaunay)
//    {
//        Foam::fvMeshWithCentres mesh
//        (
//            Foam::IOobject
//            (
//                Foam::fvMesh::defaultRegion,
//                runTime.timeName(),
//                runTime,
//                Foam::IOobject::MUST_READ
//            )
//        );

//        // Identify the patch to dualise and then extrude
//        const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);

//        // Create the dual points and faces
//        pointField dualPoints(ePatch.faceCentres());
//        faceList   dualFaces(ePatch.nPoints());
//        const labelListList& pointFaces = ePatch.pointFaces();
//            
//        forAll(dualFaces, ip)
//        {
//            const labelList& f = pointFaces[ip];
//            dualFaces[ip].setSize(f.size());
//            forAll(f, i)
//            {
//                dualFaces[ip][i] = f[i];
//            }
//            
//            // Change order if necessary
//            labelList& df = dualFaces[ip];
//            for(bool noSwaps = false; !noSwaps;)
//            {
//                const Foam::point& a = dualPoints[df[0]];
//                noSwaps = true;
//                for(label i = 1; i < df.size()-1; i++)
//                {
//                    const Foam::point& b = dualPoints[df[i]];
//                    const Foam::point& c = dualPoints[df[i+1]];            
//                    // Flip direction if necessary
//                    if ((((b-a)^(c-b)) & (a+b+c)) < 0) 
//                    {
//                        noSwaps = false;
//                        label dfi = df[i];
//                        df[i] = df[i+1];
//                        df[i+1] = dfi;
//                    }
//                }
//            }
//        }
//    
//        newPatch = PrimitivePatch<face, List, pointField>
//        (
//            dualFaces, dualPoints
//        );
//    }
    
    // Read mesh extrusion properties from earthProperties
    IOdictionary earthProperties
    (
        IOobject
        (
            "earthProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    extrudeModels::linearRadial radialExtrude(earthProperties);
    extrudedMesh newMesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime
        ),
        newPatch,
        radialExtrude
    );
    newMesh.removeBoundary();
    DynamicList<polyPatch*> newPatches(2);
    newPatches.append
    (
//        new polyPatch
        new emptyPolyPatch
        (
            "originalPatch", newPatch.size(),
            newMesh.nInternalFaces(), 0, newMesh.boundaryMesh()
        )
    );
    newPatches.append
    (
        new emptyPolyPatch
//        new polyPatch
        (
            "otherSide", newPatch.size(), 
            newMesh.nInternalFaces() + newPatch.size(), 1, newMesh.boundaryMesh()
        )
    );
    newMesh.addPatches(newPatches);
    newMesh.write();
    
    // Writing the cell and face centres of the primal mesh
    if (!writeDelaunay)
    {
        const scalar r1 = radialExtrude.Rsurface();
        const scalar r2 = radialExtrude.Router();
        //const scalar Rmid = Foam::sqrt((sqr(r1) + r1*r2 + sqr(r2))/3.);
        const scalar Rmid = 0.5*(r1 + r2);
        pointIOField newCellCentres
        (
            IOobject("cellCentres", runTime.timeName(), "polyMesh", newMesh),
            triPoints()/mag(triPoints())*Rmid
        );
        newCellCentres.write();

        // find the cross-over points between Voronoi and Delaunay meshes
        // and override the dual face centres
        pointIOField faceCentres
        (
            IOobject("faceCentres", runTime.timeName(), "polyMesh", newMesh),
            pointField(newMesh.nFaces())
        );
        for(label ie = 0; ie < newMesh.nInternalFaces(); ie++)
        {
            // find the intersection between the points at either end of the
            // edge and the Delaunay edges
            faceCentres[ie] = newCellCentres[newMesh.faceOwner()[ie]]
                            + newCellCentres[newMesh.faceNeighbour()[ie]];
            faceCentres[ie] *= Rmid/mag(faceCentres[ie]);
        }

        label m = newMesh.nInternalFaces();
        forAll(newCellCentres, i)
        {
            faceCentres[m + i] = r1/Rmid*newCellCentres[i];
        }
        m = newMesh.nInternalFaces() + newCellCentres.size();
        forAll(newCellCentres, i)
        {
            faceCentres[m + i] = r2/Rmid*newCellCentres[i];
        }

        faceCentres.write();        
    }
}
    

void Foam::VoronoiSphereMeshing::writePatch(const fileName& fName) const
{
    Info<< "Writing patch faces to " << fName << nl << endl;
    OFstream str(fName);
    str << label(number_of_vertices()-1) << "\n(" << endl;

    Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    vit++;
    do
    {
        str << setprecision(12) << topoint(vit->point()) << nl;
        vit++;
    } while (vit != finite_vertices_end());
    str << ")\n" << endl;

    str << nl << label(number_of_finite_cells()) << nl << '(' << endl;

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        FixedList<Foam::point,3> pts;
        FixedList<label,3> fac;
        label ip = 0;

        for(label i = 0; i < 4; i++)
        {
            label cvi = cit->vertex(i)->index();
            if (cvi != 0)
            {
                if (ip >= 3)
                {
                    FatalErrorIn("VoronoiSphereMeshing::writePatch")
                    << "cell does not have exatly 3 vertices on the surface of "
                        << "the sphere.\nVertices are at locations:\n"
                        << topoint(cit->vertex(0)->point()) << " with mag "
                        << mag(topoint(cit->vertex(0)->point())) << nl
                        << topoint(cit->vertex(1)->point()) << " with mag "
                        << mag(topoint(cit->vertex(1)->point())) << nl
                        << topoint(cit->vertex(2)->point()) << " with mag "
                        << mag(topoint(cit->vertex(2)->point())) << nl
                        << topoint(cit->vertex(3)->point()) << " with mag "
                        << mag(topoint(cit->vertex(3)->point())) << nl
                        << exit(FatalError);
                }

                fac[ip] = cvi-1;
                pts[ip] = topoint(cit->vertex(i)->point());
                ip++;
            }
        }
        if ((pts[0] & (pts[1] ^ pts[2])) > 0)
        {
            str << "3(" << fac[0] << ' ' << fac[2] << ' ' << fac[1] << ")\n";
        }
        else
        {
            str << "3(" << fac[0] << ' ' << fac[1] << ' ' << fac[2] << ")\n";
        }
    }
    str << ')' << endl;
}


void Foam::VoronoiSphereMeshing::writeDualPatch(const fileName& fName) const
{
    Info<< "Writing dual patch faces to " << fName << nl << endl;
    OFstream str(fName);
    const pointField pts = dualPoints();
    str << pts;

    str << nl << label(number_of_vertices()-1) << nl << '(' << endl;

    Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    vit++;
    do
    {
        // Circulate around the facets incident to this edge (non infinite)
        Cell_circulator ccirc = incident_cells(edgeBetween(vit,vCentral()));
        Cell_circulator cStart = ccirc;

        label nCells = 0;
        OStringStream output;
        output << '(';

        //str << 1 + label(std::distance(++ccirc, ccirc)) << '(';
        ccirc = cStart;
        do
        {
            if (sphericalCell(ccirc))
            {
                label celli = ccirc->cellIndex();
                if (celli >= 0)
                {
                    nCells++;
                    output << ' ' << celli;
                }
            }
        } while (--ccirc != cStart);
        output << ')';
        str << nCells << output.str().c_str() << nl;

    } while (++vit != finite_vertices_end());
    str << ')' << endl;
}


// ************************************************************************* //
