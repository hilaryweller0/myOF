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

\*----------------------------------------------------------------------------*/

#include "VoronoiSphereMeshing.H"
#include <vector>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VoronoiSphereMeshing::doubleResolution()
{
    Info << "doubleResolution:" << endl;

    label nPointsToInsert = 0;
    label nPointsTotal = 0;

    // Lists of vertices to be inserted
    pointField pointsToInsert(label(number_of_cells()*3));

    // Loop through all the edges on the sphere and create a new vertex
    for
    (
        Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        // Find the two vertices connected by this edge
        Vertex_handle vi = eit->first->vertex(eit->second);
        Vertex_handle vj = eit->first->vertex(eit->third);

        // Only for edges on sphere
        if (vi != vCentral() && vj != vCentral())
        {
            // points at either end
            Foam::point pi = topoint(vi->point());
            Foam::point pj = topoint(vj->point());
            pointsToInsert[nPointsToInsert++] = 0.5*(pi + pj);
        }
    }

    // insert the points to be inserted
    pointsToInsert.setSize(nPointsToInsert);
    insertPoints(pointsToInsert);

    Info << " Adding  " << nPointsToInsert << " points." << endl;

    reIndex();

    nPointsTotal += nPointsToInsert;
}


Foam::label Foam::VoronoiSphereMeshing::removePoints()
{
    Info << "removePoints:" << endl;

    label nPointsRemove = 0;

    // Loop through all vertices>0 and check resolution
    Finite_vertices_iterator vih = finite_vertices_begin();
    vih++;
    do
    {
        // Current point
        Foam::point defVert0 = topoint(vih->point());
        scalar dx2Ideal = sqr(initialPoints.requiredResolution(defVert0));
        scalar dxMean2 = 0;
        label np = 0;

        // All vertices connected to current vertex
        std::list<Vertex_handle> incidentVerts;
        incident_vertices(vih, std::back_inserter(incidentVerts));
        // loop around the incident vertices to calculate distances
        for
        (
            std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
            vin != incidentVerts.end();
            ++vin
        )
        {
            // only for vertices on the sphere
            if (*vin != infinite_vertex() && *vin != vCentral())
            {
                Foam::point otherPoint = topoint((*vin)->point());
                scalar dx2 = magSqr(defVert0 - otherPoint);
                dxMean2 += dx2;
                np++;
            }
        }
        //if (dxMean2 < np*0.7*dx2Ideal)
        //if (3*dxMean2 < np*dx2Ideal)
        if (dxMean2 < 0)//i.e. turned off
        {
            Finite_vertices_iterator vir = vih;
            vih--;
            remove(vir);
            nPointsRemove++;
        }
    } while(++vih != finite_vertices_end());

    Info << " Removed " << nPointsRemove << " points" << endl;

    if (nPointsRemove > 0) reIndex();

    if (number_of_vertices() <= 4)
    {
        WarningIn("VoronoiSphereMeshing::removePoints")
            << "only " << label(number_of_vertices())
            << " points left" << endl;
    }

    return nPointsRemove;
}


Foam::label Foam::VoronoiSphereMeshing::addPoints()
{
    Info << "addPoints:" << endl;

    label nPointsToInsert = 0;
    label nPointsTotal = 0;
    scalar maxDistance = 0; //distance of an edge from cente of refined region
    Info<< "Reading setResolutionDict" << endl;
    IOdictionary resDict
    (
        IOobject
        (
            "setResolutionDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    const scalar lonCentre(readScalar(resDict.lookup("lonCentre"))*pi/180.);
    const scalar latCentre(readScalar(resDict.lookup("latCentre"))*pi/180.);
    
    // the radius of the sphere
    const scalar a = mag(mesh.C()[0]);
    
    // the Cartesian form of the centre of the refined region
    const point rCentre
    (
        a*Foam::cos(lonCentre)*Foam::cos(latCentre),
        a*Foam::sin(lonCentre)*Foam::cos(latCentre),
        a*Foam::sin(latCentre)
    );
    
    
    do
    {
        nPointsToInsert = 0;

        // Lists of vertices to be inserted and cells to insert them near
        pointField pointsToInsert(label(number_of_cells()*3));
        List<Vertex_handle> insertNearVerts(label(number_of_cells()*3));

        // Loop through all the edges to find the furthest refining point
        for
        (
            Finite_edges_iterator eit = finite_edges_begin();
            eit != finite_edges_end();
            ++eit
        )
        {   // Find the two vertices on this edge
            Vertex_handle vi = eit->first->vertex(eit->second);
            Vertex_handle vj = eit->first->vertex(eit->third);
            
            //check only edges on sphere
            if (vi != vCentral() && vj != vCentral())
            {   Info << eit;
                
                //convert into FOAM points
                Foam::point pi = topoint(vi->point());
                Foam::point pj = topoint(vj->point());

                //midpoint of edge
                Foam::point midPoint = 0.5*(pi + pj);
                
                //dx2 and dx2Ideal for edge
                scalar dx2 = magSqr(pi - pj);
                scalar dx2Ideal = sqr(initialPoints.requiredResolution(midPoint));
                
                //If resolution too coarse, record distance from centre of refinement, if biggest yet
                if (dx2 > 16.0/9.0*dx2Ideal && mag(midPoint-rCentre) > maxDistance)
                {
                    maxDistance = mag(midPoint-rCentre)   
                }
            }
        }
        
        // Loop through all facets and mark for refinement those edges that belong to a face < MaxDistance from rCentre
        for
        (
            Finite_facets_iterator fit = finite_facets_begin();
            fit != finite_facets_end();
            ++fit
        )
        {   Vertex_handle vk = fit->first->vertex(fit->second);
            //check only facets on sphere,
            if (vk == vCentral())
            {  
                //convert into FOAM points
                Foam::point pi = topoint(vi->point());
                Foam::point pj = topoint(vj->point());

                //midpoint of edge
                Foam::point midPoint = 0.5*(pi + pj);
                
                //dx2 and dx2Ideal for edge
                scalar dx2 = magSqr(pi - pj);
                scalar dx2Ideal = sqr(initialPoints.requiredResolution(midPoint));
                
                //If resolution too coarse, record distance from centre of refinement, if biggest yet
                if (dx2 > 16.0/9.0*dx2Ideal && magSqr(midPoint-rCentre) > maxDistance)
                {
                    maxDistance = magSqr(midPoint-rCentre)   
                }
            }
        }
        
        // insert the points to be inserted
        pointsToInsert.setSize(nPointsToInsert);
        insertPoints(pointsToInsert, insertNearVerts);

        Info << " Adding  " << nPointsToInsert << " points." << endl;

        if (nPointsToInsert > 0) reIndex();

        nPointsTotal += nPointsToInsert;
    } while (nPointsToInsert > 0);

    return nPointsTotal;
}


Foam::label Foam::VoronoiSphereMeshing::removeBadShapes()
{
    Info << "removeBadShapes: ";

    label nPointsRemove = 0;
//    label nPointsToInsert = 0;

    // Lists of vertices to be removed (and inserted)
//    pointField pointsToInsert(label(number_of_cells()*3/2));
//    List<Vertex_handle> insertNearVerts(label(number_of_cells()*3/2));
    boolList removePoint(number_of_vertices(), false);

    // Loop through all the vertices to remove those with 4 neighbours and
    // for those with >=8 neighbours, remove the closest vertex

    Finite_vertices_iterator vih = finite_vertices_begin();
    vih++;
    do
    {
        // All finite vertices connected to current vertex
        std::list<Vertex_handle> incidentVerts;
        incident_vertices(vih, std::back_inserter(incidentVerts));
        label nInc = label
        (
            std::distance(incidentVerts.begin(),incidentVerts.end())
        )-2;

        if (nInc <= 4)
        {
            removePoint[vih->index()] = true;
            nPointsRemove++;
        }
        else if (nInc >= 8)
        {
            // Current point
            Foam::point defVert0 = topoint(vih->point());
            // Find the shortest edge
            std::list<Vertex_handle>::iterator vinclose = incidentVerts.begin();
            scalar delta2 = magSqr(defVert0);

            for
            (
                std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
                vin != incidentVerts.end();
                ++vin
            )
            {
                if (*vin != infinite_vertex() && *vin != vCentral())
                {
                    scalar dist2 = magSqr(defVert0 - topoint((*vin)->point()));
                    if (dist2 < delta2)
                    {
                        vinclose = vin;
                        delta2 = dist2;
                    }
                }
            }
            removePoint[(*vinclose)->index()] = true;
            nPointsRemove++;
        }
    } while(++vih != finite_vertices_end());

    // remove the points to be removed
    nPointsRemove = removePoints(removePoint);

//    // insert the points to be inserted
//    pointsToInsert.setSize(nPointsToInsert);
//    insertPoints(pointsToInsert, insertNearVerts); //, false);

    if (nPointsRemove > 0)// || nPointsToInsert > 0)
    {
        Info << "Removed " << nPointsRemove << " points" << endl;
            //<< " , inserted "<< nPointsToInsert << " points" << endl;
        reIndex();
    }
    else Info << "No bad shapes" << endl;

    return nPointsRemove;// + nPointsToInsert;
}


void Foam::VoronoiSphereMeshing::enlargeAngles(const label maxEnlargeAnglesIters)
{
    label nMoved = -1;

    for(label iter = 0; iter < maxEnlargeAnglesIters && nMoved!=0; iter++)
    {
        nMoved = 0;
        bool vertexMoved = false;
        
        // loop through all the points and move them if necessary
        Finite_vertices_iterator vih = finite_vertices_begin();
        for(vih++; vih != finite_vertices_end(); vih++)
        {
            // Repeatedly find the largest angle>=90 at this vertex
            // and move the vertex to the mid-point of the opposite edge
            vertexMoved = false;
            bool vertexOk = false;
            for (label itv = 0; itv < maxEnlargeAnglesIters && !vertexOk; itv++)
            {
                Foam::point defVert0 = topoint(vih->point());
                Foam::point vertNew = defVert0;
                scalar minCosAngle = 1.;
                
                // circulate around the vertices to find the maximum angle
                // All vertices connected to current vertex
                std::list<Vertex_handle> incidentVerts;
                incidentVertices(vih, incidentVerts);
                
                Foam::point p0 = topoint(incidentVerts.back()->point());
                Foam::point p1 = p0;
                for
                (
                    std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
                    vin != incidentVerts.end();
                    vin++
                )
                {
                    p1 = topoint((*vin)->point());
                    scalar cosAngle = ((p1 - defVert0)&(p0 - defVert0))
                                     /(mag(p1 - defVert0)*mag(p0 - defVert0));
                    if (cosAngle <= 0 && cosAngle < minCosAngle)
                    {
                        vertNew = 0.5*(p0 + p1);
                        minCosAngle = cosAngle;
                    }
                    p0 = p1;
                }
                
                // Move the vertex if necessary
                vertexOk = (minCosAngle > 0);
                if (!vertexOk)
                {
                    vertexMoved = true;
                    movePoint(vih, toPoint(unitVector(vertNew)*radius));
                }
            }
            if (vertexMoved) nMoved++;
        }
        Info << "enlargeAngles iteration " << iter << " moved " << nMoved
             << " points" << endl;
    }
}


void Foam::VoronoiSphereMeshing::TomitaSprings
(
    const scalar relaxation
)
{
    // Move all the points at the same time so store the new locations
    pointField newPosition(number_of_vertices());

    // Loop through all vertices>0 to calculate the vertex displacements
    label ip = 0;
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        // Current point
        Foam::point defVert0 = topoint(vih->point());
        
        // Spring force to accumulate
        vector force(0.,0.,0.);
        
        // All vertices connected to current vertex
        std::list<Vertex_handle> incidentVerts;
        incident_vertices(vih, std::back_inserter(incidentVerts));

        // loop around the incident vertices to calculate the force
        for
        (
            std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
            vin != incidentVerts.end();
            ++vin
        )
        {
            // only for vertices on the sphere
            if
            (
                *vin != infinite_vertex() && *vin != vCentral()
            )
            {
                Foam::point otherPoint = topoint((*vin)->point());
                vector dx = otherPoint - defVert0;

                // relaxed spring length
                scalar dxi = initialPoints.requiredResolution
                (
                    0.5*(defVert0 + otherPoint)
                );

                force += unitVector(dx)*(mag(dx) - beta*dxi);
            }
        }

        newPosition[ip++] = radius*unitVector(defVert0 + relaxation*force);
    }

    // Move all the points at the same time
    ip = 0;
    vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        movePoint(vih, toPoint(newPosition[ip++]));
    }    

    // Calculate the global spring energy for writing out
    scalar totalEnergy = 0;
    
    // Loop through all the edges on the sphere
    for
    (
        Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        // Find the two vertices connected by this edge and their indices
        const Vertex_handle vi = eit->first->vertex(eit->second);
        const Vertex_handle vj = eit->first->vertex(eit->third);
        
        // Only operate on these vertices if they are on the sphere
        if (vi != vCentral() && vj != vCentral())
        {
            // The points at these vertices
            const Foam::point pi = topoint(vi->point());
            const Foam::point pj = topoint(vj->point());
            // midpoint of the edge (also the crossover point)
            const vector dx = pi - pj;
            
            scalar dxi = initialPoints.requiredResolution
            (
                0.5*(pi + pj)
            );
            totalEnergy += sqr(mag(dx) - beta*dxi);
        }
    }
    
    Info << "Total spring energy = " << totalEnergy << endl;
}

void Foam::VoronoiSphereMeshing::LaplacianSmooth()
{
    // Move all the points at the same time so store the new locations
    pointField newPosition(number_of_vertices(), Foam::point(0.,0.,0.));

    // Initialise the RMS displacement
    scalar rmsDisp = 0;

    // Loop through all vertices>0 to calculate the vertex displacements
    label ip = 0;
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        // Current point
        Foam::point defVert0 = topoint(vih->point());
        
        // All vertices connected to current vertex
        std::list<Vertex_handle> incidentVerts;
        incident_vertices(vih, std::back_inserter(incidentVerts));

        // loop around the incident vertices to calculate the new position
        label nv = 0;
        for
        (
            std::list<Vertex_handle>::iterator vin = incidentVerts.begin();
            vin != incidentVerts.end();
            ++vin
        )
        {
            // only for vertices on the sphere
            if
            (
                *vin != infinite_vertex() && *vin != vCentral()
            )
            {
                nv++;
                Foam::point otherPoint = topoint((*vin)->point());
                newPosition[ip] += otherPoint;
            }
        }
        newPosition[ip] = unitVector(newPosition[ip]/nv)*radius;
        
        rmsDisp += magSqr(defVert0 - newPosition[ip]);
        ip++;
    }

    // Move all the points at the same time
    ip = 0;
    vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        movePoint(vih, toPoint(newPosition[ip++]));
    }

    Info << "RMS displacement = " << sqrt(rmsDisp) << endl;
}


void Foam::VoronoiSphereMeshing::LloydIteration()
{
    // Move all the vertices at the same time so store the new positions
    pointField newPositions(number_of_vertices());

    // Initialise the RMS displacement
    scalar rmsDisp = 0;

    // Loop through all vertices>0
    label ip = 0;
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        // Current point
        Foam::point defVert0 = topoint(vih->point());
        
        // Voronoi cell centre to move the point to
        Foam::point vCentre(0.,0.,0.);
        label nVpoints = 0;
        // remember the previous point of the triangle
        Foam::point prevPt(0,0,0);
        
        // accumulate the area of the Voronoi cell and the centres of the
        // sub-triangles including the mesh density
        scalar sumArho = 0;
        vector sumAcrho(0,0,0);
        
        // First circulate around adjacent points to find approximate centre
        
        // Circulate around the cells incident to this edge (non infinite)
        Cell_circulator ccirc = incident_cells(edgeBetween(vih,vCentral()));
        const Cell_circulator cStart = ccirc;
        do
        {
            nVpoints++;
            prevPt = circumcenter(ccirc);
            vCentre = vCentre + prevPt;
        } while (++ccirc != cStart);
        vCentre /= scalar(nVpoints);
        
        // Circulate around again to find the area weighted Voronoi cell centre
        ccirc = cStart;
        do
        {
            Foam::point currPt = circumcenter(ccirc);
            
            //Point c = CGAL::centroid(currPt, prevPt, vCentre);
            Foam::point c = (currPt + prevPt + vCentre)/3.;
            scalar a = mag((currPt - prevPt) ^ (vCentre - prevPt));
            scalar rho = 1/pow(initialPoints.requiredResolution(c),4);
            
            sumArho += a*rho;
            sumAcrho = sumAcrho + a*rho*c;
            
            prevPt = currPt;
        } while (++ccirc != cStart);
        
        vCentre = radius*unitVector(sumAcrho/(sumArho + VSMALL));

        // Store the new Voronoi cell centre
        newPositions[ip++] = vCentre;
        //movePoint(vih, vCentre);

        // Sum the displacement
        rmsDisp += magSqr(defVert0 - vCentre);
    }
    
    // Move all the points at the same time
    ip = 0;
    vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        movePoint(vih, toPoint(newPositions[ip++]));
    }    

    Info << "RMS displacement = " << sqrt(rmsDisp) << endl;
}


void Foam::VoronoiSphereMeshing::HRify()
{
    // Make the mesh more Heikes Randall. ie move edge cross-over points to
    // mid-way between Voronoi vertices. This is done on the Delaunay grid
    // by moving the the problem Delaunay edge towards the further point, 
    // along the existing edges

    // Move all the vertices at the same time so store the displacements
    pointField disp(number_of_vertices(), Foam::point::zero);

    // Initialise the RMS displacement
    scalar rmsDisp = 0;

    // Loop through all the edges on the sphere
    for
    (
        Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        // Find the two vertices connected by this edge and their indices
        const Vertex_handle vi = eit->first->vertex(eit->second);
        const Vertex_handle vj = eit->first->vertex(eit->third);
        const label ii = vi->index()-1;
        const label ij = vj->index()-1;
        
        // Only operate on these vertices if they are on the sphere
        if (vi != vCentral() && vj != vCentral())
        {
            // The points at these vertices
            const Foam::point pi = topoint(vi->point());
            const Foam::point pj = topoint(vj->point());
            // midpoint of the edge (also the crossover point)
            const Foam::point pijmid = 0.5*(pi + pj);

            // The opposite vertices of the triangles
            Pair<Foam::point> pts = pointsOpposite(*eit);
            //Pair<label> vts = verticesOpposite(*eit);
            
            // Find which is the close vertex and which is far
            label cls = 0;
            label far = 1;
            if (magSqr(pts[cls] - pijmid) > magSqr(pts[far] - pijmid))
            {
                cls = 1;
                far = 0;
            }
            
            // The circumcentres of the triangles
            Foam::Pair<Foam::point> circum;
            circum[cls] = topoint
            (
                construct_circumcenter(toPoint(pi),toPoint(pj),toPoint(pts[cls]))
            );
            circum[far] = topoint
            (
                construct_circumcenter(toPoint(pi),toPoint(pj),toPoint(pts[far]))
            );

            // required mesh density at circumcentres
            Foam::Pair<scalar> rho;
            rho[cls] = 1/sqr(initialPoints.requiredResolution(circum[cls]));
            rho[far] = 1/sqr(initialPoints.requiredResolution(circum[far]));

            // weighted midpoint of the circumcentres
            const Foam::point cm = (rho[cls]*circum[cls] + rho[far]*circum[far])
                                    /(rho[cls] + rho[far]);

            // the displacement for pi and pj
            vector d = 0.4*(cm - pijmid);
            
            // unit vector to pull pj to pfar
            vector pjfar = unitVector(pts[far] - pj);
            
            // unit vector to pull pi to pfar
            vector pifar = unitVector(pts[far] - pi);
            
            // increment the displacements for pi and pj
            disp[ii] += (pifar & d)*pifar;
            disp[ij] += (pjfar & d)*pjfar;
        }
    }
    
    // Move all the points at the same time
    label ip = 0;
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        Foam::point newPoint = radius
                            *unitVector(topoint(vih->point()) + disp[ip]);
        rmsDisp += magSqr(topoint(vih->point()) - newPoint);
        movePoint(vih, toPoint(newPoint));
        ip++;
    }    

    Info << "RMS displacement = " << sqrt(rmsDisp) << endl;
}


Foam::label Foam::VoronoiSphereMeshing::PittewayIteration(const scalar relax)
{
    Info << "PittewayIteration(relaxation = " << relax << ") ";

    label nNonPittEdges = 0;

    // Move all the vertices at the same time so store the displacements
    pointField disp(number_of_vertices(), Foam::point::zero);
    boolList moveVertex(number_of_vertices(), false);

    // Loop through all the edges on the sphere
    for
    (
        Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        // Find the two vertices connected by this edge and their indices
        const Vertex_handle vi = eit->first->vertex(eit->second);
        const Vertex_handle vj = eit->first->vertex(eit->third);
        const label ii = vi->index()-1;
        const label ij = vj->index()-1;
        
        // Only operate on these vertices if they are on the sphere
        if (vi != vCentral() && vj != vCentral())
        {
            // The points at these vertices
            const Foam::point pi = topoint(vi->point());
            const Foam::point pj = topoint(vj->point());

            // The opposite vertices of the triangles
            Pair<Foam::point> pts = pointsOpposite(*eit);
            //Pair<label> vts = verticesOpposite(*eit);
            
            // Find if either vertex makes an angle of >=90
            label far = -1;
            scalar cosAngle = (pi - pts[0]) & (pj - pts[0])
                               /(mag(pi - pts[0])*mag(pj - pts[0]));
            if (cosAngle <= 0) far = 0;
            else
            {
                cosAngle = (pi - pts[1]) & (pj - pts[1])
                        /(mag(pi - pts[1])*mag(pj - pts[1]));
                if (cosAngle <= 0) far = 1;
            }

            // If neither is far, do nothing. Otherwise move the edge towards far
            if (far != -1)
            {
                nNonPittEdges++;
                moveVertex[ii] = true;
                moveVertex[ij] = true;
            
                // midpoint of the edge
                const Foam::point pijmid = 0.5*(pi + pj);
                // the mid-point of the circumcentres
                const Foam::point cm = 0.5*
                (
                    topoint(construct_circumcenter
                    (
                        toPoint(pi),toPoint(pj),toPoint(pts[0])
                    ))
                  + topoint(construct_circumcenter
                    (
                        toPoint(pi),toPoint(pj),toPoint(pts[1])
                    ))
                );
                // the displacement for pi and pj
                vector d = 0.5*(cm - pijmid);
                
                // increment the displacements for pi and pj
                disp[ii] += d;
                disp[ij] += d;

//                // unit vectors to pull pj and pi to pfar
//                vector pjfar = unitVector(pts[far] - pj);
//                vector pifar = unitVector(pts[far] - pi);
//                
//                // increment the displacements for pi and pj
//                disp[ii] += (pifar & d)*pifar;
//                disp[ij] += (pjfar & d)*pjfar;
            }
        }
    }
    
    // Move all the points at the same time
    label ip = 0;
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        if (moveVertex[ip])
        {
            Foam::point newPoint = radius
                           *unitVector(topoint(vih->point()) + relax*disp[ip]);
            movePoint(vih, toPoint(newPoint));
        }
        ip++;
    }    
    
    Info << " modified " << nNonPittEdges << " non-Pitteway edges" << endl;

    return nNonPittEdges;
}


void Foam::VoronoiSphereMeshing::spatialSortRenumber()
{
    Info << "Renumbering" << endl;    
    
    // Put all the points on a std::vector for sorting
    std::vector<Point> VoronoiPoints;
    VoronoiPoints.reserve(number_of_vertices()-1);
    
    Finite_vertices_iterator vih = finite_vertices_begin();
    for(vih++; vih != finite_vertices_end(); vih++)
    {
        VoronoiPoints.push_back(vih->point());
    }
    
    // sort the points
//    CGAL::hilbert_sort
//    (
//        VoronoiPoints.begin(),
//        VoronoiPoints.end(),
//        CGAL::Hilbert_sort_median_policy()
//    );
    CGAL::spatial_sort
    (
        VoronoiPoints.begin(),
        VoronoiPoints.end()
    );
    
    // remove points with old order from the triangulation and insert new ones
    removePoints(boolList(number_of_vertices(), true));
    //insert(VoronoiPoints.begin(), VoronoiPoints.end());
    label nVert = number_of_vertices();
    // Add the points and index them
    for
    (
        std::vector<Point>::iterator vpit = VoronoiPoints.begin();
        vpit != VoronoiPoints.end();
        vpit++
    )
    {
        insert(*vpit)->index() = nVert++;
    }
    
    // Renumbering the triangular cells base on the new vertex numbering
    // First give the cells numbering based on the iterator
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

    // Store the new cell indices before they are assigned
    labelList newCellIndices(number_of_finite_cells(), -1);
    
    // Dual faces for finding cell numbers from vertices
    const faceList vertexCells = dualFaces();
    
    // New cell indices
    label ic = 0;
    
    // Loop through the VoronoiPoints and find the newCellIndices
    std::vector<Point>::iterator vpit = VoronoiPoints.begin();
    for(label ip = 0; ip < label(number_of_vertices())-1; ip++)
    {
        // number the cells incident to this vertex if not already done
        for(label ivc = 0; ivc < vertexCells[ip].size(); ivc++)
        {
            label oldCell = vertexCells[ip][ivc];
            if (newCellIndices[oldCell] == -1)
            {
                newCellIndices[oldCell] = ic++;
            }
        }
        vpit++;
    }
    
    // Renumber the cells
    dualVerti = 0;
    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        cit++
    )
    {
        if (sphericalCell(cit))
        {
            cit->cellIndex() = newCellIndices[dualVerti++];
        }
    }
}

// ************************************************************************* //
