{
// Local value of pv
const volScalarField pv1 = linear ? 
                           volScalarField((TRiSK::curl(vS) + f)/H) :
                           pv;

const dimensionedScalar l1h = fvc::domainIntegrate(mag(h - hInit))/hI1sum;
const dimensionedScalar l2h = Foam::sqrt
(
    fvc::domainIntegrate(sqr(h - hInit))/hI2sum
);
const dimensionedScalar mass = fvc::domainIntegrate(h);
const dimensionedScalar lih = max(mag(h - hInit))/hInitMax;

const dimensionedScalar l1pv = fvc::domainIntegrate(mag(pv1 - pvInit))/pvI1sum;
const dimensionedScalar l2pv = Foam::sqrt
(
    fvc::domainIntegrate(sqr(pv1 -pvInit))/pvI2sum
);
const dimensionedScalar lipv = max(mag(pv1 - pvInit))/pvInitMax;

const dimensionedScalar l1U = TRiSK::domainIntegrate(mag(Uf - UfInit))/UfI1sum;
const dimensionedScalar l2U = Foam::sqrt
(
    TRiSK::domainIntegrate((Uf - UfInit) & (Uf - UfInit))/UfI2sum
);
const dimensionedScalar liU = max(mag(Uf - UfInit))/UfInitMax;

energyOld = energy;
energy = linear? H*(TRiSK::KE(TRiSK::circToFlux(vS), vS) + magg*(0.5*h + h0))
               : h*KE + magg*h*(0.5*h + h0);
dimensionedScalar energySum = fvc::domainIntegrate(energy)/Vtot;

surfaceScalarField CoriolisEnergy = phi*dualMesh.signedDualMap(pvFlux)
                                   /(mesh.depthf()*mesh.faceVol());

volScalarField divhuDiff = fvc::div(phi) - h*fvc::div(phi/hf)
               - 2*TRiSK::faceToCellMap(phi/(hf*mesh.magSf())*fvc::snGrad(h));
dimensionedScalar sumDivhuDiff = fvc::domainIntegrate(divhuDiff)/Vtot;

const dimensionedScalar PVmean = fvc::domainIntegrate(pv1)/Vtot;
const dimensionedScalar enstrophy = linear ?
    0.5*fvc::domainIntegrate(sqr(pv1)*H)/Vtot
  : 0.5*fvc::domainIntegrate(sqr(pv1)*hv)/Vtot;

os << runTime.timeName() << "  " << l2h.value() << "  " << lih.value()
   << "  " << ((energySum - energyInit)/energyInit).value()
   << "  " << (enstrophy - enstrophyInit).value()/max(mag(enstrophyInit).value(), SMALL)
   << "  " << (PVmean - PVinit).value()
//   << "  " << KEdoubleTime.value()/(3600*24*365)
   << "  " << ((mass - massInit)/massInit).value()
   << "  " << (TRiSK::domainIntegrate(CoriolisEnergy)/Vtot).value()
   << "  " << sumDivhuDiff.value()
   << "  " << l1h.value()
   << "  " << l1pv.value()
   << "  " << l2pv.value()
   << "  " << lipv.value()
   << "  " << l1U.value()
   << "  " << l2U.value()
   << "  " << liU.value()
   << endl;

//const dimensionedScalar dualMass = fvc::domainIntegrate(hv);
//const dimensionedScalar faceMass = TRiSK::domainIntegrate(hf);
//const dimensionedScalar faceMassDual = TRiSK::domainIntegrate(mesh.dualMap(hf));
//Info << "Dual mass error = " << ((dualMass - mass)/mass).value() << nl
//     << "Face mass error = " << ((faceMass - mass)/mass).value() << nl
//     << "Face mass error = " << ((faceMassDual - mass)/mass).value() << endl;
}

