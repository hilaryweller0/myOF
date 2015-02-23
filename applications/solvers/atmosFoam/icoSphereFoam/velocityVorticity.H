// update full velocity (for diagnostics and for CLUST)
Uf = TRiSK::divReconstruct3d(uS);

// Calculation of PV flux and KE

// vertical componenet of pv
pv = f;
if (!linear) pv += TRiSK::curl(mesh.dualFluxMap(vS));
pvzf = dualMesh.dualMap(fvc::interpolate(pv));

//pvzf = ff;
//if (!linear)
//{
//    pvzf += dualMesh.dualMap
//    (
//        fvc::interpolate(TRiSK::curl(mesh.dualFluxMap(vS)))
//    );
//}

// 3d pv
pvf = 2*Omega;
if (!linear) pvf += TRiSK::curl3d(vS);

// replace vertical part of pv
//pvh = pvf - (pvf & mesh.rHatf())*mesh.rHatf();
//pvf = pvh + pvzf*mesh.rHatf();

// only use vertical part of pv
//pvf = pvzf*djperp;

// Calculate pv flux
pvFlux = (pvf ^ Uf) & mesh.Sf();
//if (conserveEnergy)
//{
//    pvFlux = 0.5*
//    (
//        pvFlux
//      + dualMesh.dualFluxMap(TRiSK::perp(pvzf*vS))
////      + ((pvh ^ Uf) & mesh.Sf())
//    );
//}

// Kinetic energy
KE = TRiSK::KE(uS, vS);
