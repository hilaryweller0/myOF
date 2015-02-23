Info << "Initialising diagnostics file errorDiags.dat\n" << endl;

fileName invs = args.rootPath() / args.caseName() / "errorDiags.dat";
OFstream os(invs);
os << "#time   l2h   lInfh      normalisedEnergy   normalisedEnstrophy    PVdiff  normlisedMass CoriolisEnergy PGenergy  l1h l1pv l2pv lInfpv l1U l2U lInfU"
   << endl;
//curlGradh zonalAsym latMaxh"

// Initialising diagnostics
volScalarField energy = linear 
                  ? H*(TRiSK::KE(TRiSK::circToFlux(vS), vS) + magg*(0.5*h + h0))
                  : h*TRiSK::KE(TRiSK::circToFlux(vS), vS) + magg*h*(0.5*h +h0);
volScalarField energyOld = energy;

const volScalarField hInit = h;
const volScalarField pvInit = linear ? volScalarField((TRiSK::curl(vS) + f)/H) : pv;
const surfaceVectorField UfInit = Uf;

const dimensionedScalar Vtot = sum(mesh.V());
const dimensionedScalar hmin("hmin", dimLength, SMALL);
const dimensionedScalar massInit = max(fvc::domainIntegrate(h), hmin*min(mesh.V()));
const dimensionedScalar hI1sum = max(fvc::domainIntegrate(mag(hInit)), hmin*min(mesh.V()));
const dimensionedScalar hI2sum = max(fvc::domainIntegrate(sqr(hInit)), sqr(hmin)*min(mesh.V()));
const dimensionedScalar hInitMax = max(max(mag(hInit)), hmin);
const dimensionedScalar PVinit = fvc::domainIntegrate(pvInit)/Vtot;

const dimensionedScalar pvI1sum = fvc::domainIntegrate(mag(pvInit));
const dimensionedScalar pvI2sum = fvc::domainIntegrate(sqr(pvInit));
const dimensionedScalar pvInitMax = max(mag(pvInit));

const dimensionedScalar umin("umin", dimVelocity, SMALL);
const dimensionedScalar UfI1sum = max(TRiSK::domainIntegrate(mag(UfInit)), umin*min(mesh.V()));
const dimensionedScalar UfI2sum = max(TRiSK::domainIntegrate(UfInit&UfInit), umin*umin*min(mesh.V()));
const dimensionedScalar UfInitMax = max(max(mag(UfInit)), umin);

const dimensionedScalar energyInit = fvc::domainIntegrate(energy)/Vtot;
//const dimensionedScalar KEInit = fvc::domainIntegrate(KE)/Vtot;
//dimensionedScalar KEdoubleTime("KEdoubleTime", dimTime, scalar(0));
const dimensionedScalar enstrophyInit = linear ?
    0.5*fvc::domainIntegrate(sqr(pvInit)*H)/Vtot
  : 0.5*fvc::domainIntegrate(sqr(pvInit)*hv)/Vtot;

//const dimensionedScalar dualVol = sum(dualMesh.V());
//const dimensionedScalar faceVol = sum(mesh.faceVol());
//const dimensionedScalar faceVolDual = sum(dualMesh.faceVol());
//Info << "Dual volume error = " << ((dualVol - Vtot)/Vtot).value() << nl
//     << "Face volume error = " << ((faceVol - Vtot)/Vtot).value() << nl
//     << "Dual face volume error = " << ((faceVolDual - Vtot)/Vtot).value()
//     << endl;

