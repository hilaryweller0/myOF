namespace Foam
{
// Moves the dual patch points to reduce skewness
void optimiseDualPoints
(
    pointField& oldDualPoints,
    const faceList& dualFaces, 
    const polyPatch& primalPatch,
    const label nIters
)
{
    PrimitivePatch<face, List, pointField> dualPatch(dualFaces, oldDualPoints);

    const pointField& dualPoints = dualPatch.localPoints();
    const pointField& primalPoints = primalPatch.localPoints();
    
    // First find map from dual edges to primal edges
    labelList edgeMap(dualPatch.nEdges(), -1);
    for(label ie = 0; ie < dualPatch.nEdges(); ie++)
    {
        const edge& dualEdge = dualPatch.edges()[ie];
        // points either end of this edge (same numbers as faces of the primal)
        const label id0 = dualEdge[0];
        const label id1 = dualEdge[1];
        const label ip0 = dualPatch.meshPoints()[id0];
        const label ip1 = dualPatch.meshPoints()[id1];
        
        for
        (
            label iei = 0;
            iei < primalPatch.faceEdges()[ip0].size()
         && edgeMap[ie] == -1;
            iei++
        )
        {
            label edgeI = primalPatch.faceEdges()[ip0][iei];
            // loop around edges of primal face id1
            for
            (
                label iej = 0;
                iej < primalPatch.faceEdges()[ip1].size()
             && edgeMap[ie] == -1;
                iej++
            )
            {
                if (primalPatch.faceEdges()[ip1][iej] == edgeI)
                {
                    edgeMap[ie] = edgeI;
                }
            }
        }
        // Check that an edge has been found
        if (edgeMap[ie] == -1)
        {
            FatalErrorIn("optimiseDualPoints")
                << "cannot find a primal edge which crosses over dual edge "
                << ie << abort(FatalError);
        }
    }

    // Declare new points to move the dual points to
    pointField newPoints(dualPoints);
    // displacements for the points
    vectorField disp(newPoints.size(), vector::zero);
    // Radius of the dual patch (to stay the same)
    const scalar radius = mag(newPoints[0]);
    
    // Move the points nIters times
    for(label it = 0; it < nIters; it++)
    {
        // Loop over all the points of the dual
        for(label id = 0; id < dualPatch.nPoints(); id++)
        {
            // Initialise displacement for the dual point
            disp[id] = vector::zero;
            
            // Loop over all the edges of the point and increment the displacement
            for(label ie = 0; ie < dualPatch.pointEdges()[id].size(); ie++)
            {
                // primal and dual edges
                label ide = dualPatch.pointEdges()[id][ie];
                const edge& dEdge = dualPatch.edges()[ide];
                const edge& pEdge = primalPatch.edges()[edgeMap[ide]];
               
                // primal and dual points
                const point& dp0 = newPoints[dEdge[0]];
                const point& dp1 = newPoints[dEdge[1]];
                const point& pp0 = primalPoints[pEdge[0]];
                const point& pp1 = primalPoints[pEdge[1]];
                
                // primal edge centre (with the correct magnitude)
                const point edgeCentre = radius*unitVector(pp0 + pp1);
                
                // primal and dual edge intersection
                plane pPlane(pp0 ^ pp1);
                plane dPlane(dp0 ^ dp1);
                plane::ray r = pPlane.planeIntersect(dPlane);
                point intersect = radius*unitVector(r.dir())
                                *sign(r.dir()&edgeCentre);
                
                // displacement is from intersection to primal centre
                disp[id] += edgeCentre - intersect;
            }
            disp[id] /= scalar(dualPatch.pointEdges()[id].size());
        }
        newPoints += 0.2*disp;
    }
    
    // rearrnage the newPoints into the global point order
    forAll(oldDualPoints, ip)
    {
        oldDualPoints[ip] = newPoints[dualPatch.meshPointMap()[ip]];
    }
}
}
