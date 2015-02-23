Info << "Initialising diagnostics file errorDiags.dat\n" << endl;

fileName invs = args.rootPath() / args.caseName() / "errorDiags.dat";
OFstream os(invs);
os << "#time   l2h   lInfh      normalisedEnergy   normalisedEnstrophy    PVdiff  normlisedMass CoriolisEnergy PGenergy  l1h l1pv l2pv lInfpv l1U l2U lInfU"
   << endl;
//curlGradh zonalAsym latMaxh"

// Initialising diagnostics
volScalarField energy = linear ? H*(TRiSK::ke(un) + magg*(0.5*h + h0))
                               : h*TRiSK::ke(un) + magg*h*(0.5*h + h0);
volScalarField energyOld = energy;

const volScalarField hInit = h;
const volScalarField pvInit = linear ? volScalarField((TRiSK::curl(un) + f)/H) : pv;
const surfaceVectorField UfInit = Uf;

const dimensionedScalar Vtot = sum(mesh.V());
const dimensionedScalar massInit = fvc::domainIntegrate(h);
const dimensionedScalar hI1sum = fvc::domainIntegrate(mag(hInit));
const dimensionedScalar hI2sum = fvc::domainIntegrate(sqr(hInit));
const dimensionedScalar hInitMax = max(mag(hInit));
const dimensionedScalar PVinit = fvc::domainIntegrate(pvInit)/Vtot;

const dimensionedScalar pvI1sum = fvc::domainIntegrate(mag(pvInit));
const dimensionedScalar pvI2sum = fvc::domainIntegrate(sqr(pvInit));
const dimensionedScalar pvInitMax = max(mag(pvInit));

const dimensionedScalar UfI1sum = TRiSK::domainIntegrate(mag(UfInit));
const dimensionedScalar UfI2sum = TRiSK::domainIntegrate(UfInit&UfInit);
const dimensionedScalar UfInitMax = max(mag(UfInit));
Info <<  "UfI1Sum=" << UfI1sum << endl;
Info <<  "UfI2Sum=" << UfI2sum << endl;

const dimensionedScalar energyInit = fvc::domainIntegrate(energy)/Vtot;
//const dimensionedScalar KEInit = fvc::domainIntegrate(KE)/Vtot;
//dimensionedScalar KEdoubleTime("KEdoubleTime", dimTime, scalar(0));
const dimensionedScalar enstrophyInit = linear ?
    0.5*fvc::domainIntegrate(sqr(pvInit)*H)/Vtot
  : 0.5*fvc::domainIntegrate(sqr(pvInit)*TRiSK::primalToDualCellMap(h))/Vtot;

