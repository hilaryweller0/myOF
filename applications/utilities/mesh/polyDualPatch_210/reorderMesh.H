#include "fvMesh.H"
#include "polyTopoChange.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "bandCompression.H"
#include "faceSet.H"
#include "SortableList.H"

// cellOrder: old cell for every new cell
// faceOrder: old face for every new face. Ordering of boundary faces not
// changed.
autoPtr<mapPolyMesh> reorderMesh
(
    polyMesh& mesh,
    const labelList& cellOrder,
    const labelList& faceOrder
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));

    faceList newFaces(reorder(reverseFaceOrder, mesh.faces()));
    labelList newOwner
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceOwner())
        )
    );
    labelList newNeighbour
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceNeighbour())
        )
    );

//    // Check if any faces need swapping.
//    forAll(newNeighbour, faceI)
//    {
//        label own = newOwner[faceI];
//        label nei = newNeighbour[faceI];

//        if (nei < own)
//        {
//            newFaces[faceI] = newFaces[faceI].reverseFace();
//            Swap(newOwner[faceI], newNeighbour[faceI]);
//            Info << "Swapping owner " << own << " with neighbour " << nei
//                 << " of face " << faceI << endl;
//        }
//    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());

    forAll(patches, patchI)
    {
        patchSizes[patchI] = patches[patchI].size();
        patchStarts[patchI] = patches[patchI].start();
        oldPatchNMeshPoints[patchI] = patches[patchI].nPoints();
        patchPointMap[patchI] = identity(patches[patchI].nPoints());
    }

    mesh.resetPrimitives
    (
        Xfer<pointField>::null(),
        xferMove(newFaces),
        xferMove(newOwner),
        xferMove(newNeighbour),
        patchSizes,
        patchStarts,
        true
    );

    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,                       //const polyMesh& mesh,
            mesh.nPoints(),             // nOldPoints,
            mesh.nFaces(),              // nOldFaces,
            mesh.nCells(),              // nOldCells,
            identity(mesh.nPoints()),   // pointMap,
            List<objectMap>(0),         // pointsFromPoints,
            faceOrder,                  // faceMap,
            List<objectMap>(0),         // facesFromPoints,
            List<objectMap>(0),         // facesFromEdges,
            List<objectMap>(0),         // facesFromFaces,
            cellOrder,                  // cellMap,
            List<objectMap>(0),         // cellsFromPoints,
            List<objectMap>(0),         // cellsFromEdges,
            List<objectMap>(0),         // cellsFromFaces,
            List<objectMap>(0),         // cellsFromCells,
            identity(mesh.nPoints()),   // reversePointMap,
            reverseFaceOrder,           // reverseFaceMap,
            reverseCellOrder,           // reverseCellMap,
            labelHashSet(0),            // flipFaceFlux,
            patchPointMap,              // patchPointMap,
            labelListList(0),           // pointZoneMap,
            labelListList(0),           // faceZonePointMap,
            labelListList(0),           // faceZoneFaceMap,
            labelListList(0),           // cellZoneMap,
            pointField(0),              // preMotionPoints,
            patchStarts,                // oldPatchStarts,
            oldPatchNMeshPoints         // oldPatchNMeshPoints
        )
    );
}

