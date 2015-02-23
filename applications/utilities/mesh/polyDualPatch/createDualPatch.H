// Identify the patch to dualise and then extrude
const polyPatch& ePatch(mesh.boundaryMesh()["originalPatch"]);

// Create the dual points and faces
pointField dualPoints(ePatch.faceCentres());
//pointField dualPoints(mesh.cellCentres());
faceList   dualFaces(ePatch.nPoints());
const labelListList& pointFaces = ePatch.pointFaces();

forAll(dualFaces, ip)
{
    const labelList& f = pointFaces[ip];
    dualFaces[ip].setSize(f.size());
    forAll(f, i)
    {
        dualFaces[ip][i] = f[i];
    }
    
    // Change order if necessary
    labelList& df = dualFaces[ip];
    for(bool noSwaps = false; !noSwaps;)
    {
        const point& a = dualPoints[df[0]];
        noSwaps = true;
        for(label i = 1; i < df.size()-1; i++)
        {
            const point& b = dualPoints[df[i]];
            const point& c = dualPoints[df[i+1]];
            // Flip direction if necessary
            if ((((b-a)^(c-b)) & (a+b+c)) < 0) 
            {
                noSwaps = false;
                label dfi = df[i];
                df[i] = df[i+1];
                df[i+1] = dfi;
            }
        }
    }
}

