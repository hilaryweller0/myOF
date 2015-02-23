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
#include "extrudedMesh.H"
#include "linearRadial.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField Foam::VoronoiSphereMeshing::dualPoints
(
    const faceList& faces
) const
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

    // Collapse points which are close together (if needed)
    label nCollapse = 0;
    if (mag(minEdgeLength) > SMALL)
    {
        forAll(faces, faci)
        {
            const face f = faces[faci];
            // first circulate around the face and find the longest edge
            label iprev = f.size()-1;
            scalar maxEdgeLength = 0;
            for(label ip = 0; ip < f.size(); ip++)
            {
                if (mag(dPoints[f[ip]] - dPoints[f[iprev]]) > maxEdgeLength)
                {
                    maxEdgeLength = mag(dPoints[f[ip]] - dPoints[f[iprev]]);
                }
                iprev = ip;
            }
            
            // circulate around the face again and check edge lengths
            iprev = f.size()-1;
            for(label ip = 0; ip < f.size(); ip++)
            {
                if
                (
                    mag(dPoints[f[ip]] - dPoints[f[iprev]])
                         < minEdgeLength*maxEdgeLength
                )
                {
                    Foam::point midPoint = radius*unitVector
                    (
                        dPoints[f[ip]] + dPoints[f[iprev]]
                    );
                    dPoints[f[ip]] = midPoint;
                    dPoints[f[iprev]] = midPoint;
                    nCollapse ++;
                }
                iprev = ip;
            }
        }
    }
    if (nCollapse > 0)
    {
        Info << "Collapsed " << nCollapse << " pairs of points" << endl;
    }

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
    Info << "Writing the primal mesh with centres" << endl;
    
    // Create the primal mesh patch to extrude
    const faceList dualFace = dualFaces();
    PrimitivePatch<face, List, pointField> newPatch
    (
        dualFace, dualPoints(dualFace)
    );

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
    
    // remove the empty zeroth boundary
    newMesh.removeBoundary();
    DynamicList<polyPatch*> newPatches(2);
    if(radialExtrude.nLayers() == 1)
    {
        newPatches.append
        (
            new emptyPolyPatch
            (
                "originalPatch", newPatch.size(),
                newMesh.nInternalFaces(), 0, newMesh.boundaryMesh(),
                "empty"
            )
        );
        newPatches.append
        (
            new emptyPolyPatch
            (
                "otherSide", newPatch.size(), 
                newMesh.nInternalFaces()+newPatch.size(),1,newMesh.boundaryMesh(),
                "empty"
            )
        );
    }
    else
    {
        newPatches.append
        (
            new polyPatch
            (
                "originalPatch", newPatch.size(),
                newMesh.nInternalFaces(), 0, newMesh.boundaryMesh(),
                "patch"
            )
        );
        newPatches.append
        (
            new polyPatch
            (
                "otherSide", newPatch.size(), 
                newMesh.nInternalFaces()+newPatch.size(),1,newMesh.boundaryMesh(),
                "patch"
            )
        );
    }
    newMesh.addPatches(newPatches);
    
    newMesh.write();
    
    const pointField triPts = triPoints();
    
    // Writing the cell and face centres
    scalarList radiiC(radialExtrude.nLayers());
    scalarList radiiCf(radialExtrude.nLayers()+1);
    for(label k = 0; k < radialExtrude.nLayers()+1; k++)
    {
        radiiCf[k] = mag(radialExtrude
        (
            newPatch.points()[0], newPatch.pointNormals()[0], k
        ));
    }
    for(label k = 0; k < radialExtrude.nLayers(); k++)
    {
        radiiC[k] = 0.5*(radiiCf[k] + radiiCf[k+1]);
    }
    
//    const scalar r1 = radialExtrude.Rsurface();
//    const scalar r2 = radialExtrude.Router();
//    //const scalar Rmid = Foam::sqrt((sqr(r1) + r1*r2 + sqr(r2))/3.);
//    const scalar Rmid = 0.5*(r1 + r2);
    
    pointIOField newCellCentres
    (
        IOobject("cellCentres", runTime.timeName(), "polyMesh", newMesh),
        newMesh.nCells()
    );
    label cellI = 0;
    // if use Voronoi generation points as cell centres:
    if (VoronoiCentres)
    {
        for (label layer = 0; layer < radiiC.size(); layer++)
        {
            for(label i = 0; i < triPts.size(); i++)
            {
                newCellCentres[cellI++] = triPts[i]
                                      /mag(triPts[i])*radiiC[layer];
            }
        }
    }
    else // use centroids
    {
        for (label layer = 0; layer < radiiC.size(); layer++)
        {
            for(label i = 0; i < triPts.size(); i++)
            {
                newCellCentres[cellI++] = newPatch.faceCentres()[i]
                                 /mag(newPatch.faceCentres()[i])*radiiC[layer];
            }
        }
    }
    
    newCellCentres.write();

    // find the cross-over points between Voronoi and Delaunay meshes
    // and override the dual face centres
        pointIOField faceCentres
        (
            IOobject("faceCentres", runTime.timeName(), "polyMesh", newMesh),
            pointField(newMesh.nFaces())
        );

    const label nEdges = newPatch.nEdges();
    const label nCols = triPts.size();

    // First the vertical faces (all internal)
    for(label ie = 0; ie < nEdges; ie++)
    {
        // find the intersection between the points at either end of the
        // edge and the Delaunay edges
        Foam::point fc = newCellCentres[newMesh.faceOwner()[ie]]
                       + newCellCentres[newMesh.faceNeighbour()[ie]];
        fc /= mag(fc);
        for(label k = 0; k < radiiC.size(); k++)
        {
            label faceI = ie + k*(nEdges + nCols);
            faceCentres[faceI] = fc*radiiC[k];
        }
    }
    
    // Next the horizontal faces
    for(label ic = 0; ic < triPts.size(); ic++)
    {
        Foam::point fc = triPts[ic]/mag(triPts[ic]);
        
        // first the internal horizontal faces
        for(label k = 0; k < radiiCf.size()-2; k++)
        {
            label faceI = ic + k*nCols + (k+1)*nEdges;
            faceCentres[faceI] = fc*radiiCf[k+1];
        }
        
        // Bottom and top faces
        faceCentres[ic+newMesh.nInternalFaces()] = fc*radiiCf[0];
        faceCentres[ic+nCols+newMesh.nInternalFaces()] = fc*radiiCf.last();
    }

    faceCentres.write();
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
    const pointField pts = dualPoints(dualFaces());
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
