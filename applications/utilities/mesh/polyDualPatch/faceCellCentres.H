    // Writing the cell and face centres (of the dual)
    scalarList radiiC(radialExtrude.nLayers());
    scalarList radiiCf(radialExtrude.nLayers()+1);
    for(label k = 0; k < radialExtrude.nLayers()+1; k++)
    {
        radiiCf[k] = mag(radialExtrude
        (
            dualPatch.points()[0], dualPatch.pointNormals()[0], k
        ));
    }
    for(label k = 0; k < radialExtrude.nLayers(); k++)
    {
        radiiC[k] = 0.5*(radiiCf[k] + radiiCf[k+1]);
    }

    Info << "Writing the dual mesh cell centres\n" << endl;
    
    pointIOField newCellCentres
    (
        IOobject("cellCentres", runTime.timeName(), "polyMesh", dualMesh),
        dualMesh.nCells()
    );
    pointField patchFaceCentres = ePatch.localPoints()/mag(ePatch.localPoints());
    label cellI = 0;
    for (label layer = 0; layer < radiiC.size(); layer++)
    {
        for(label i = 0; i < nCols; i++)
        {
            newCellCentres[cellI++] = patchFaceCentres[i]*radiiC[layer];
        }
    }
    newCellCentres.write();

    Info << "Calculating the primal and dual mesh face centres\n" << endl;
    
    pointIOField newFaceCentres
    (
        IOobject
        (
            "faceCentres",
            runTime.findInstance(mesh.meshDir(), "points"),
            "polyMesh",
            mesh
        ),
        mesh.faceCentres()
    );
    
    pointIOField dualFaceCentres
    (
        IOobject("faceCentres", runTime.timeName(), "polyMesh", dualMesh),
        pointField(dualMesh.nFaces())
    );
    
//    if (Cftype == CROSS)
//    {
//        Info << "Using primal-dual cross-over points" << endl;
        // First the vertical faces (all internal) - dual-primal cross-over points
        for(label ie = 0; ie < nEdges; ie++)
        {
            for(label k = 0; k < radiiC.size(); k++)
            {
                label faced = ie + k*(nEdges + nCols);
                label faceI = dualPrimalFaceMap[faced];
                //newFaceCentres[faced] = mesh.faceCentres()[faceI];
                const point& Cpo = mesh.C()[mesh.owner()[faceI]];
                const point& Cpn = mesh.C()[mesh.neighbour()[faceI]];
                const point& Cdo = newCellCentres[dualMesh.faceOwner()[faced]];
                const point& Cdn = newCellCentres[dualMesh.faceNeighbour()[faced]];
                
                if (mag((Cpo^Cpn) ^ (Cdo^Cdn)) < SMALL)
                {
                    WarningIn("faceCellCentres")
                         << "primal and dual edges are parallel" << nl;
                    Info << "Cpo = " << Cpo << nl
                         << "Cpn = " << Cpn << nl
                         << "Cdo = " << Cdo << nl
                         << "Cdn = " << Cdn << nl;
                }
                
                plane pFace(Cpo^Cpn);
                plane dFace(Cdo^Cdn);

                plane::ray r = pFace.planeIntersect(dFace);
                newFaceCentres[faceI] = unitVector(r.dir())*sign(r.dir() & Cpo)
                                       *radiiC[k];
                dualFaceCentres[faced] = newFaceCentres[faceI];
                const point & Cf = newFaceCentres[faceI];
                
                // Check that primal and dual acrually cross
                if (sign((Cf - Cdo) & (Cf - Cdn)) > 0)
                {
                    WarningIn("faceCellCenrres.H") << "primal face " << faceI
                    << " does not overlap with dual face " << faced << endl;
                }
            }
        }
//    }
//    else if (Cftype == PRIMAL)
//    {
//        Info << "Using the primal face centres" << endl;
//        for(label ie = 0; ie < nEdges; ie++)
//        {
//            for(label k = 0; k < radiiC.size(); k++)
//            {
//                label faced = ie + k*(nEdges + nCols);
//                label faceI = dualPrimalFaceMap[faced];
//                dualFaceCentres[faced] = newFaceCentres[faceI];
//            }
//        }
//    }
//    else if (Cftype == DUAL || Cftype == ERROR)
//    {
//        FatalErrorIn("polyDualPatch") << "Cf type should be cross (for "
//            << "cross-over) or primal for primal face centres. Dual face "
//            << "centres not yet implemented" << exit(FatalError);
//    }

    // Next the horizontal faces
    for(label ic = 0; ic < nCols; ic++)
    {
        Foam::point fc = patchFaceCentres[ic];
        
        // first the internal horizontal faces
        for(label k = 0; k < radiiCf.size()-2; k++)
        {
            label faced = ic + k*nCols + (k+1)*nEdges;
            dualFaceCentres[faced] = fc*radiiCf[k+1];
        }
        
        // Bottom and top faces
        dualFaceCentres[ic+dualMesh.nInternalFaces()] = fc*radiiCf[0];
        dualFaceCentres[ic+nCols+dualMesh.nInternalFaces()] =fc*radiiCf.last();
    }

    newFaceCentres.write();
    dualFaceCentres.write();
    
    // Check distance between dual and primal faces
    scalarList faceDist(dualMesh.nInternalFaces(), scalar(0));
    forAll(faceDist, faced)
    {
        if (dualPrimalFaceMap[faced] != -1)
        {
            faceDist[faced] = mag
            (
                newFaceCentres[dualPrimalFaceMap[faced]]
              - dualFaceCentres[faced]
            );
        }
    }
    Info << "Maximum distance between primal and dual face centres = "
         << max(faceDist)<< endl;
    

