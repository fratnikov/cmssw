#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

using namespace std;
using namespace reco;

namespace {
  enum HGCType {HGCEE=0,HGCHEF=1,HGCHEB=2,HGCUNDEF=3};

  HGCType hitType (const PFRecHit& recHit) {
    DetId id = recHit.detId();
    if (id.det() == DetId::Forward) {
      if (id.subdetId() == ForwardSubdetector::HGCEE) return HGCType::HGCEE;
      if (id.subdetId() == ForwardSubdetector::HGCHEF) return HGCType::HGCHEF;
      if (id.subdetId() == ForwardSubdetector::HGCHEB) return HGCType::HGCHEB;
    }
    return HGCUNDEF;
  }

  int HGCRecHitLayer (const PFRecHit& recHit) {
    HGCType type = hitType (recHit);
    if (type == HGCEE) {
      //      cout << "HGCRecHitLayer-> " << HGCEEDetId(recHit.detId()) << endl;
      return HGCEEDetId (recHit.detId()).layer();
    }
    if (type == HGCHEF || type == HGCHEB) {
      //      cout << "HGCRecHitLayer-> " << HGCHEDetId(recHit.detId()) << endl;
      return HGCHEDetId (recHit.detId()).layer();
    }
    return -1; // undefined
  }

  class RecalibrationConstants {
  public:
    RecalibrationConstants (const edm::ParameterSet& pSet) 
      : mipValueInGeV (pSet.getParameter<double>("MipValueInGeV")),
	coef_a (pSet.getParameter<double>("effMip_to_InverseGeV_a")),	
	coef_b (pSet.getParameter<double>("effMip_to_InverseGeV_b")),	
	coef_c (pSet.getParameter<double>("effMip_to_InverseGeV_c")),
	weights (pSet.getParameter<vector<double> >("weights"))
    { 
    }

    double calibrateRecHitEnergy (const PFRecHit& recHit) const {
      double eta = recHit.position().eta();
      double energy_MIP = recHit.energy()/mipValueInGeV;
      double eCorr = weights[HGCRecHitLayer(recHit)]*energy_MIP;//std::cosh(eta);
      double effMIP_to_InvGeV = coef_a/(1.0 + std::exp(-coef_c - coef_b*std::cosh(eta)));
      // cout << "calibrateRecHitEnergy-> layer "<<HGCRecHitLayer(recHit)
      // 	   << " weight:a:b:c " << weights[HGCRecHitLayer(recHit)]<<':'<<coef_a<<':'<<coef_b<<':'<<coef_c
      // 	   << " energy_MIP:eCorr:effMIP_to_InvGeV:result " <<energy_MIP<<':'<<eCorr<<':'<<effMIP_to_InvGeV<<':'<<eCorr/effMIP_to_InvGeV
      // 	   << endl;
      return eCorr/effMIP_to_InvGeV;

    }
  private:
    double mipValueInGeV;
    double coef_a;
    double coef_b;
    double coef_c;
    vector<double> weights;
    
  };

}

class PFRecHitCalibrator: public  edm::EDProducer {
public:
  PFRecHitCalibrator (const edm::ParameterSet&);
  virtual ~PFRecHitCalibrator ();
  virtual void produce (edm::Event&, const edm::EventSetup&) override;
private:
  edm::InputTag fInputLabel;
  const RecalibrationConstants constantsHGC [3];
};

PFRecHitCalibrator::PFRecHitCalibrator (const edm::ParameterSet& iConfig)
  : fInputLabel (iConfig.getParameter<edm::InputTag>("input")),
    constantsHGC {iConfig.getParameterSet("HGCEEElectronEnergyCalibrator"),
                  iConfig.getParameterSet("HGCHEFHadronicEnergyCalibrator"),
                  iConfig.getParameterSet("HGCHEBHadronicEnergyCalibrator")
    }
{
  produces<RecoPFRecHitRefCandidateCollection>();
}

PFRecHitCalibrator::~PFRecHitCalibrator () {
}

void PFRecHitCalibrator::produce (edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  edm::Handle<RecoPFRecHitRefCandidateCollection> inputHits;
  iEvent.getByLabel(fInputLabel, inputHits);
  //  std::cout << "PFRecHitCalibrator::produce-> input " << inputHits->size() << std::endl;
  std::auto_ptr<RecoPFRecHitRefCandidateCollection> calibratedHits  (new RecoPFRecHitRefCandidateCollection ());
  for (auto inHit : *inputHits) {
    RecoPFRecHitRefCandidate outHit (inHit);
    const PFRecHit& recHit = *(inHit.pfRecHit());
    HGCType type = hitType (recHit);
    if (type != HGCUNDEF) {
      int index = int (type);
      double newEnergy = constantsHGC[index].calibrateRecHitEnergy (recHit);
      outHit.setEnergy (newEnergy);
    }
    if (1 || fabs(inHit.eta())<0.1) {
      // cout << "PFRecHitCalibrator-> det "<<int(DetId(recHit.detId()).det())<<':'<<DetId(recHit.detId()).subdetId()<<" eta:phi:energy " << inHit.eta()<<':'<<inHit.phi()<<':'<<inHit.energy()
      // 	  << " -> " <<outHit.energy()
      // 	  << endl;
    }
    calibratedHits->push_back (outHit);
  }
  iEvent.put(calibratedHits);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFRecHitCalibrator);
