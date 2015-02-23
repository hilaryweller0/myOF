// Control Parameters for VoronoiSphereMeshing

//- Write out the mesh every iteration
const Switch writeIntermediate;

//- Write the Voronoi generation points as the cell centres or the centroids
const Switch VoronoiCentres;

//- Relaxation factor for the iteration
const scalar relaxationFactorStart;

//- Relaxation factor for the final iteration
const scalar relaxationFactorEnd;

//- Relaxation factor for iteration number it (out of maxIt)
scalar relaxationFactor(const label it, const label maxIt)
{
    return relaxationFactorStart
        + it/max(maxIt,scalar(1))*(relaxationFactorEnd - relaxationFactorStart);
}

// Double the global resolution before anything else
const Switch doubleGlobalResolution;

// Add or points to achieve the required resolution?
const Switch addRemovePoints;

// Change the topology or keep fixed
const Switch changeTopology;

// Should we attempt to fix bad shapes
const Switch fixBadShapes;

// Number of fake Lloyd iterations
const label nFakeLloyds;

// Number of Lloyd iterations
const label nLloydIterations;

// Number of Laplacian smoothing iterations
const label nLaplacianSmooths;

// Number of iterations to make the mesh Heikes and Randall
const label nHRiterations;

// Maximum number of iterations to make the mesh Pitteway
const label maxPittewayIterations;

// Number of Tomita spring relaxations
const label nTomitaSprings;

//- Ratio between natural spring length and required resolution
const scalar beta;

//- smooth all points or just irregular ones
const Switch smoothAll;

//- renumber all the points using a CGAL spatial sort
const Switch spatialSort;

//- Collapse edges in the Voronoi dual less than minEdgeLength* the largest edge
//- of that polygonal cell
const scalar minEdgeLength;


