namespace Foam
{
void removeTriangles(faceList& faces, pointField& points)
{
    PrimitivePatch<face, List, pointField> patch(faces, points);
    
    const pointField& lPoints = patch.localPoints();
    const faceList& lFaces = patch.localFaces();
    
    // Find points to remove and the points to replace them with
    labelList replacePointWith(lPoints.size(), -1);
    DynamicList<point> newPoints;
    boolList keepFace(lFaces.size(), true);
    
    // Loop through all the faces and find triangles
    forAll(lFaces, faci)
    {
        const face& f = lFaces[faci];
        if (f.size() == 3)
        {
            keepFace[faci] = false;
            // Find the triangle centre
            point fCtr = 1./3.*(lPoints[f[0]] + lPoints[f[1]] + lPoints[f[2]]);
            newPoints.append(fCtr);
            // mark all the triange points to be removed
            replacePointWith[f[0]] = newPoints.size()-1;
            replacePointWith[f[1]] = newPoints.size()-1;
            replacePointWith[f[2]] = newPoints.size()-1;
        }
    }
    
    // Create new point list with additional points
    pointField allPoints(lPoints);
    allPoints.append(newPoints);
    
    // Assign new faceList
    faceList newFaces(lFaces.size()-newPoints.size());
    label newFaceI = 0;
    for(label faci = 0; faci < lFaces.size(); faci++)
    {
        if (keepFace[faci])
        {
            // Check if any of the points are to go
            bool keepWholeFace = true;
            for(label fi = 0; fi < lFaces[faci].size() && keepWholeFace; fi++)
            {
                if(replacePointWith[lFaces[faci][fi]] != -1)
                {
                    keepWholeFace = false;
                }
            }
            
            if (keepWholeFace) newFaces[newFaceI] = lFaces[faci];
            else
            {
                bool replaced = false;
                newFaces[newFaceI].setSize(lFaces[faci].size());
                label nfi = 0;
                for(label fi = 0; fi < lFaces[faci].size(); fi++)
                {
                    if(replacePointWith[lFaces[faci][fi]] == -1)
                    {
                        newFaces[newFaceI][nfi] = lFaces[faci][fi];
                        nfi++;
                    }
                    else if (!replaced)
                    {
                        newFaces[newFaceI][nfi]
                           = lPoints.size() + replacePointWith[lFaces[faci][fi]];
                        nfi++;
                        replaced = true;
                    }
                }
                if (nfi != lFaces[faci].size()) newFaces[newFaceI].setSize(nfi);
            }
            newFaceI++;
        }
    }
    
    // replace the original faces and points
    faces = newFaces;
    points = allPoints;
}
}
