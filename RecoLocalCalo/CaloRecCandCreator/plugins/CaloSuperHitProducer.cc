#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "DataFormats/RecoCandidate/interface/CaloSuperHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

using namespace std;
using namespace reco;

namespace {
  const double TWOPI = 2.*3.1415927; 

  ostream& operator<< (ostream& out, DetId det) {
    return out<<int(det.det())<<':'<<det.subdetId();
  }

}

class CaloSuperHitProducer: public  edm::EDProducer {
public:
  CaloSuperHitProducer (const edm::ParameterSet&);
  virtual ~CaloSuperHitProducer ();
  virtual void produce (edm::Event&, const edm::EventSetup&) override;
  virtual void beginRun (const edm::Run&, const edm::EventSetup&) override;
  size_t etaIndex (double eta) const;
  size_t phiIndex (double eta, double phi) const;
  pair<double,double> binEtaDEta (size_t iEta) const;
  pair<double,double> binPhiDPhi (size_t iEta, size_t iPhi) const;
private:
  size_t localEta (double eta) const;
  edm::InputTag fInputLabel;
  vector<unsigned> fPhiBins;
  vector<double> fEtaBins;
};

CaloSuperHitProducer::CaloSuperHitProducer (const edm::ParameterSet& iConfig)
  : fInputLabel (iConfig.getParameter<edm::InputTag>("input")),
    fPhiBins (iConfig.getParameter<vector<unsigned> >("phiBins")),
    fEtaBins (iConfig.getParameter<vector<double> >("etaBins"))
{
  if (fPhiBins.size() != fEtaBins.size()) {
    std::cout << "CaloSuperHitProducer::CaloSuperHitProducer-> mismatch between Eta and Phi sizes!" << endl;
    fPhiBins.resize (fEtaBins.size(), fPhiBins.back());
  }
  produces<vector<CaloSuperHit> >();
}

CaloSuperHitProducer::~CaloSuperHitProducer () {
}

void CaloSuperHitProducer::beginRun (const edm::Run& iRun, const edm::EventSetup& iSetup) {
}

void CaloSuperHitProducer::produce (edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  edm::Handle<vector<CaloRecHitCandidate> > inputHits;
  iEvent.getByLabel(fInputLabel, inputHits);
  std::cout<<"CaloSuperHitProducer::produce-> " << inputHits->size()<<" input hits"<<endl;
  std::auto_ptr<vector<CaloSuperHit> > cleanedSuperHits  (new vector<CaloSuperHit> ());
  // 1 loop fill energies
  map < pair<unsigned,unsigned>, CaloSuperHit > superHits;
  for (auto& inHit : *inputHits) {
    double eta = inHit.eta();
    double phi = inHit.phi();
    // fill SuperHits
    pair<size_t, size_t> key (etaIndex (eta), phiIndex (eta, phi));
    CaloSuperHit& superHit = superHits[key];
    if (superHit.etaPhiSquare () <= 0) {
      pair<double,double> etaDEta = binEtaDEta (key.first);
      pair<double,double> phiDPhi = binPhiDPhi (key.first, key.second);
      superHit = CaloSuperHit(etaDEta.first, phiDPhi.first, etaDEta.second, phiDPhi.second);
    }
    superHit.addHit (inHit.detid(), inHit.p4());
  }
  vector<double> allEnergies (fEtaBins.size(), 0.);
  vector<double> allEnergies2 (fEtaBins.size(), 0.);
  for (const auto& sh : superHits) {
    const CaloSuperHit& superHit = sh.second;
    size_t iEta = localEta(superHit.etaCenter ());
    double p = superHit.p();
    allEnergies[iEta] += p;
    allEnergies2[iEta] += p*p;
  }
  // 3 loop remove energy from SuperHits
  for (auto& sh : superHits) {
    CaloSuperHit& superHit = sh.second;
    double pNew = superHit.p();
    size_t iEta = localEta(superHit.etaCenter ());
    double nCells = TWOPI/superHit.dPhi();
    double meanP = allEnergies[iEta]/nCells;
    double meanP2 = allEnergies2[iEta]/nCells;
    double sigma = sqrt (meanP2-meanP*meanP);
    pNew -= (meanP+sigma);
    // if (iEta == 15) {
    //   cout<<"CaloSuperHitProducer::produce-> ieta:ieta:iphi "<<iEta<<':'<<sh.first.first<<':'<<sh.first.second<<" nCells:"<<nCells
    // 	  <<" eta:phi "<<superHit.etaCenter()<<':'<<superHit.phiCenter()/TWOPI*360
    // 	  <<" orig/meanP/sqrt(meanP2)/meanP2/sigma/cleaned: "<<superHit.p()<<'/'<<meanP<<':'<<sqrt(meanP2)<<'/'<<meanP2<<'>'<<sigma<<':'<<pNew<<endl;
    // }
    if (pNew > 0) {
      double scale = pNew / superHit.p();
      math::PtEtaPhiMLorentzVector newP4 (scale*superHit.pt(),superHit.eta(),superHit.phi(),0);
      superHit.setP4 (newP4);
      cleanedSuperHits->push_back (superHit);
    }
  }
  iEvent.put(cleanedSuperHits);
}

size_t CaloSuperHitProducer::localEta (double eta) const {
  auto ptr = lower_bound (fEtaBins.begin(), fEtaBins.end(), fabs(eta));
  if (ptr == fEtaBins.end()) return fEtaBins.size()-1;
  return ptr - fEtaBins.begin();
}

size_t CaloSuperHitProducer::etaIndex (double eta) const {
  size_t iEta = localEta (fabs (eta));
  size_t nBins = fEtaBins.size();
  if (eta >= 0) return iEta+nBins;
  else return nBins-iEta-1;
}

size_t CaloSuperHitProducer::phiIndex (double eta, double phi) const {
  size_t nPhiBins = fPhiBins [localEta (fabs (eta))];
  if (phi < 0) phi += TWOPI;
  return size_t (floor (phi/TWOPI*nPhiBins));
}

pair<double,double> CaloSuperHitProducer::binEtaDEta (size_t iEta) const {
  pair<double,double> result;
  size_t nBins = fEtaBins.size();
  size_t localEta = (iEta >= nBins) ? iEta - nBins : nBins - iEta - 1; 
  if (localEta == 0) { 
    result.first = result.second = 0.5*fEtaBins[0];
  }
  else {
    result.first = 0.5*(fEtaBins[localEta] + fEtaBins[localEta-1]);
    result.second = 0.5*(fEtaBins[localEta] - fEtaBins[localEta-1]);
  }
  return result;
}

pair<double,double> CaloSuperHitProducer::binPhiDPhi (size_t iEta, size_t iPhi) const {
  size_t nBins = fEtaBins.size();
  size_t localEta = iEta >= nBins ? iEta - nBins : nBins - iEta - 1;
  double dPhi = TWOPI/fPhiBins[localEta];
  return pair<double,double> (dPhi*(iPhi+0.5), 0.5*dPhi);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloSuperHitProducer);
