namespace Foam
{
// Moves the dual patch points to reduce skewness
labelList findEdgeMap
(
    PrimitivePatch<face, List, pointField>& dualPatch, 
    const polyPatch& primalPatch
)
{
    const label nEdges = dualPatch.nEdges();
    labelList edgeMap(nEdges, -1);
    for(label ie = 0; ie < nEdges; ie++)
    {
        // faces of the dual
        const label i = dualPatch.edgeFaces()[ie][0];
        const label j = dualPatch.edgeFaces()[ie][1];

        // edges of the primal associated with these primal faces
        const labelList& ipEdges = primalPatch.pointEdges()[i];
        const labelList& jpEdges = primalPatch.pointEdges()[j];

        // In the primal patch, find the edge between points i and j
        for(label iei = 0; iei < ipEdges.size() && edgeMap[ie] == -1; iei++)
        {
            for(label jej = 0; jej < jpEdges.size() &&edgeMap[ie]==-1;jej++)
            {
                if (ipEdges[iei] == jpEdges[jej])
                {
                    edgeMap[ie] = ipEdges[iei];
                }
            }
        }
        
        // Check that the edge is found and is within bounds
        if (edgeMap[ie] == -1 || edgeMap[ie] >= nEdges)
        {
            FatalErrorIn("findEdgeMap")
                << "Cannot find edge in between primal patch points "
                << i << " and " << j <<  " for dual edge " << ie << nl
                << "ipEdges = " << ipEdges << nl
                << "jpEdges = " << jpEdges << nl
                << exit(FatalError);
        }
    }
    
    return edgeMap;
}
}
