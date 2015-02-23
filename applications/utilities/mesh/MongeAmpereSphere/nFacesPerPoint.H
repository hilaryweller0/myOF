// Calculate the number of faces per point (assuming that the mesh is fixed)
// for mapping gradPhi from faces to points

labelList nFacesPerPoint(mesh.nPoints(), label(0));
for(label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
{
    const face& f = mesh.faces()[faceI];
    forAll(f, fi)
    {
        nFacesPerPoint[f[fi]] ++;
    }
}

