// Additional Mesh navigation functions for VoronoiSphereMeshing

//- Edges around a vertex
//void incident_edges
//(
//    const Vertex_handle& v, std::list<Edge>& edges
//) const;

////- Edges of facets incident to a vertex which are opposite the vertex
//void opposite_edges
//(
//    const Vertex_handle& v, std::list<Edge>& edges
//) const;

//- Edge opposite to vertices v0 and v1 in cell c
Edge edgeOpposite
(
    const Cell_handle& ch,
    const Vertex_handle& v0,
    const Vertex_handle& v1
) const;

Edge edgeBetween
(
    const Vertex_handle& v0, const Vertex_handle& v1
) const;

//- The circumcenter of the face on the sphere of cell c
Foam::point circumcenter(const Cell_handle& c) const;

//- Test if the cell is a spherical cell (ie attached to the central
//  vertex)
bool sphericalCell(const Cell_handle& c) const;

//- The 2 points opposite edge e
Foam::Pair<Foam::point> pointsOpposite(const Edge& e) const;

//- The 2 vertices opposite edge e
Foam::Pair<label> verticesOpposite(const Edge& e) const;

//- The 2 finite cells neighbouring edge e and the points opposite the
//- edge of each cell
void neighbourCells
(
    Foam::Pair<Cell_handle>& neibCells,
    Foam::Pair<Foam::point>& pOut,
    const Edge& e
) const;

//- Return central vertex handle
inline Vertex_handle vCentral() const
    {return Vertex_handle(finite_vertices_begin());};

//- Set the list of incident vertices to vertex vit on the sphere
void incidentVertices
(
    const Vertex_handle& vit, std::list<Vertex_handle>& incidentVerts
) const;

