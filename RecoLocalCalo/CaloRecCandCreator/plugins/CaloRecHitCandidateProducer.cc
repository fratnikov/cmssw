//
// Producer to converts ECAL and HCAL hits into Candidates
// Author: F. Ratnikov (Maryland)
// Jan. 7, 2007
// Updated: F.Ratnikov (Nebraska)
// Dec. 3, 2014
//


#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDProducer.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

#include "RecoLocalCalo/CaloRecCandCreator/interface/CaloRecHitRecalibrator.h"

using namespace edm;
using namespace reco;
using namespace std;

namespace {

  ostream& operator<< (ostream& stream, DetId detId){
    if (detId.det() == DetId::Forward) {
      if (detId.subdetId() == ForwardSubdetector::HGCEE) {
	stream << HGCEEDetId (detId.rawId());
      }
      else {
	stream << HGCHEDetId (detId.rawId());
      }
    }
    else {
      stream<<int(detId.det())<<':'<<detId.subdetId();
    }
    return stream;
  }

  const CaloSubdetectorGeometry* getCaloSubdetectorGeometry (const edm::EventSetup & fSetup, DetId detId) {
    if (detId.det() == DetId::Forward) {
      edm::ESHandle<HGCalGeometry> hgGeometry;
      ForwardSubdetector subdet = static_cast<ForwardSubdetector>(unsigned(detId.subdetId()));
      string nameX = subdet == ForwardSubdetector::HGCEE ? "HGCalEESensitive" :
	subdet == ForwardSubdetector::HGCHEF ? "HGCalHESiliconSensitive" :
	subdet == ForwardSubdetector::HGCHEB ? "HGCalHEScintillatorSensitive" :
	"";
      fSetup.get<IdealGeometryRecord>().get(nameX, hgGeometry);
      return static_cast<const CaloSubdetectorGeometry*> (&*hgGeometry);
    }
    else {
      ESHandle<CaloGeometry> geometry;
      fSetup.get<CaloGeometryRecord>().get (geometry);
      return geometry->getSubdetectorGeometry (detId);
    }
  }

  GlobalPoint cellPosition (const CaloSubdetectorGeometry* geometry, DetId detId) {
    GlobalPoint result (0,0,0);
    if (detId.det() == DetId::Forward) {
      result = dynamic_cast <const HGCalGeometry*> (geometry)->getPosition (detId);
      //      cout<<"CellPosition-> "<<detId<<" "<<result.x()<<':'<<result.y()<<':'<<result.z()
      //	  <<" "<<result.eta()<<':'<<result.phi()<<endl;
    }
    else {
      const CaloCellGeometry* cellGeometry =geometry->getGeometry (detId);
      if (cellGeometry) {
	result = cellGeometry->getPosition ();
      }
    }
    return result;
  }
}

class CaloRecHitCandidateProducer : public edm::EDProducer {
public:
  CaloRecHitCandidateProducer( const edm::ParameterSet&); 
  virtual ~CaloRecHitCandidateProducer() {delete mRecalibrator;}
  void produce( edm::Event&, const edm::EventSetup& );
  
private:
  edm::InputTag mInputLabel;
  double mEnergyThreshold;
  CaloRecHitRecalibrator* mRecalibrator;
};

CaloRecHitCandidateProducer::CaloRecHitCandidateProducer ( const edm::ParameterSet & fConfig ) 
  :  mInputLabel (fConfig.getParameter<edm::InputTag>("input")),
     mEnergyThreshold (fConfig.getParameter<double>("energyThreshold")),
     mRecalibrator (CaloRecHitRecalibrator::newCaloRecHitRecalibrator (fConfig.getParameterSet ("recalibrator")))
{
  produces<vector<CaloRecHitCandidate> >();
}

void CaloRecHitCandidateProducer::produce( edm::Event & fEvent, const edm::EventSetup & fSetup) {
  auto_ptr<vector<CaloRecHitCandidate> > output ( new vector<CaloRecHitCandidate> );
  edm::Handle<edm::View<CaloRecHit> > caloRecHits;
  fEvent.getByLabel(mInputLabel, caloRecHits);
  DetId lastDetId (0);
  const CaloSubdetectorGeometry* geometry = 0;
  for (size_t ihit = 0; ihit < caloRecHits->size(); ++ihit) {
    const CaloRecHit& caloRecHit = (*caloRecHits)[ihit];
    DetId id = caloRecHit.detid();
    if (id.det() != lastDetId.det() || id.subdetId() != lastDetId.subdetId()) {
      geometry = getCaloSubdetectorGeometry (fSetup, id);
      lastDetId = id;
    }
    GlobalPoint pos = cellPosition (geometry, id);
    double eta = pos.eta();
    double energy = mRecalibrator->recalibrateEnergy (caloRecHit.energy(), pos.eta(), id);
    //    cout<<"CaloRecHitCandidateProducer::produce-> eta="<<eta<<':'<<pos.phi()<<" energy "<<caloRecHit.energy()<<" -> " <<energy<<endl;
    if (energy <= mEnergyThreshold) continue;
    math::RhoEtaPhiVector p( 1, eta, pos.phi());
    p *= ( energy / p.r() );
    CaloRecHitCandidate crh ( Candidate::LorentzVector( p.x(), p.y(), p.z(), energy ) );
    crh.setDetId (RefToBase<CaloRecHit>(caloRecHits, ihit));
    output->push_back( crh );
  }
  fEvent.put(output);
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloRecHitCandidateProducer);
