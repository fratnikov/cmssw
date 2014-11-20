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
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

using namespace std;
using namespace reco;

namespace {

  ostream& operator<< (ostream& out, DetId det) {
    return out<<int(det.det())<<':'<<det.subdetId();
  }

  bool sameSubdet (DetId a, DetId b) {
    return a.det() == b.det() && a.subdetId() == b.subdetId();
  }

  double xySquare (const GlobalPoint& a, const GlobalPoint& b, const GlobalPoint& c) {
    double dx = sqrt((a.x()-b.x())*(a.x()-b.x())+(a.y()-b.y())*(a.y()-b.y()));
    double dy = sqrt((c.x()-b.x())*(c.x()-b.x())+(c.y()-b.y())*(c.y()-b.y()));
    GlobalPoint center (0.5*(a.x()+c.x()), 0.5*(a.y()+c.y()), 0.5*(a.z()+c.z()));
    GlobalPoint ref1 (0.5*dx, center.perp()+0.5*dy, center.z());
    GlobalPoint ref2 (-0.5*dx, center.perp()-0.5*dy, center.z());
    double dEta = fabs (ref1.eta()-ref2.eta());
    double dPhi = fabs (ref1.phi()-ref2.phi());
    return dEta*dPhi;
  }
  
  class GenericCaloGeometry {
  public:
    GenericCaloGeometry () {};
    void addDetector (DetId det);
    void reset (const edm::EventSetup& iSetup);
    const CaloSubdetectorGeometry* getGeometry (DetId det) const;
    const vector<DetId>& getValidDetIds (DetId det) const;
    double eta (DetId det) const;
    GlobalPoint position (DetId det) const;
    double etaPhiSquare (DetId det) const;
    vector<DetId> allDetectors () const;
  private:
    map <DetId, const CaloSubdetectorGeometry*> geometries;
  };
  
  void GenericCaloGeometry::addDetector (DetId det) {
    DetId id (det.det(), det.subdetId());
    if (geometries.find (id) == geometries.end()) {
      geometries[id] = 0;
    }
  }
  
  void GenericCaloGeometry::reset (const edm::EventSetup& iSetup) {
    for (auto& pg : geometries) {
      DetId det = pg.first;
      if (det.det() == DetId::Forward) {
	string detKey;
	if (det.subdetId() == ForwardSubdetector::HGCEE) detKey = "HGCalEESensitive";
	else if (det.subdetId() == ForwardSubdetector::HGCHEF) detKey = "HGCalHESiliconSensitive";
	else if (det.subdetId() == ForwardSubdetector::HGCHEB) detKey = "HGCalHEScintillatorSensitive";
	else {
	  pg.second = 0;
	  continue;
	}
	edm::ESHandle<HGCalGeometry> hgcg;
	iSetup.get<IdealGeometryRecord>().get(detKey, hgcg);
	pg.second = static_cast<const CaloSubdetectorGeometry*> (&*hgcg);
      }
      else if (det.det() == DetId::Ecal || det.det() == DetId::Hcal) {
	edm::ESHandle<CaloGeometry> cgp;
	iSetup.get<CaloGeometryRecord> ().get (cgp);
	pg.second = cgp->getSubdetectorGeometry(det.det(), det.subdetId());
      }
    }
  }
  
  const CaloSubdetectorGeometry* GenericCaloGeometry::getGeometry (DetId det) const {
    DetId id (det.det(), det.subdetId());
    auto pg = geometries.find (id);
    if (pg != geometries.end()) return pg->second;
    return 0;
  }
  
  const vector<DetId>& GenericCaloGeometry::getValidDetIds (DetId det) const {
    const CaloSubdetectorGeometry* geometry = getGeometry(det);
    // if (det.det() == DetId::Forward) {
    //   return dynamic_cast<const HGCalGeometry*>(geometry)->getValidGeomDetIds();
    // }
    return geometry->getValidDetIds(det.det(), det.subdetId());
  }
  
  double GenericCaloGeometry::eta (DetId det) const {
    const CaloSubdetectorGeometry* geometry = getGeometry(det);
    if (geometry) {
      if (det.det() == DetId::Forward) {
	return dynamic_cast<const HGCalGeometry*>(geometry)->getPosition(det).eta();
      }
      else {
	return getGeometry(det)->getGeometry(det)->etaPos();
      }
    }
    cerr<<"GenericCaloGeometry::eta-> Uninitiated detector "<<det<<endl;
    return 0;
  }
  
  GlobalPoint GenericCaloGeometry::position (DetId det) const {
    if (det.det() == DetId::Forward) {
      return dynamic_cast<const HGCalGeometry*>(getGeometry(det))->getPosition (det);
    }
    else {
      return getGeometry(det)->getGeometry(det)->getPosition ();
    }
  }
  
  double GenericCaloGeometry::etaPhiSquare (DetId det) const {
    const CaloSubdetectorGeometry* geometry = getGeometry (det);
    if (!geometry) return 0;
    if (det.det() == DetId::Forward) {
      const HGCalGeometry* fwdGeometry = dynamic_cast<const HGCalGeometry*> (getGeometry (det));
      if (!fwdGeometry) return 0;
      vector<GlobalPoint> corners = fwdGeometry->getCorners(det);
      return xySquare (corners[0], corners[1], corners[2]);
    }
    else if (det.det() == DetId::Ecal || det.det() == DetId::Hcal) {
      return geometry->deltaPhi(det)*geometry->deltaEta(det);
    }
    return 0;
  }
  
  vector<DetId> GenericCaloGeometry::allDetectors () const {
    vector<DetId> result;
    for (auto& pg : geometries) {
      result.push_back(pg.first);
    }
    return result;
  }
  
  class LayerDetId : public DetId {
  public:
    LayerDetId (DetId id);
    LayerDetId (DetId id, unsigned int layer);
    unsigned int layer ();
  private:
    void setLayer (unsigned layer);
    void setLayer ();
  };
  
  LayerDetId::LayerDetId (DetId id) 
      : DetId (id)
    {
      setLayer ();
    }

    LayerDetId::LayerDetId (DetId id, unsigned int layer) 
      : DetId (id) 
    {
      setLayer (layer);
    }

    unsigned int LayerDetId::layer () {
      return rawId() & 0xFF;
    }

    void LayerDetId::setLayer (unsigned layer) {
      id_ = (id_ & 0xFE000000) | (layer & 0xFF);
    }

    void LayerDetId::setLayer () {
      unsigned int layer = 0;
      if (det() == Ecal) {
	if (subdetId() == EcalSubdetector::EcalBarrel) layer = 1;
	else if (subdetId() == EcalSubdetector::EcalShashlik) layer = 1;
      }
      else if (det() == Hcal) {
	if (subdetId() == HcalSubdetector::HcalBarrel || 
	    subdetId() == HcalSubdetector::HcalEndcap || 
	    subdetId() == HcalSubdetector::HcalForward){
	  layer = HcalDetId(rawId()).depth();
	}
      }
      else if (det() == Forward) {
	if (subdetId() == ForwardSubdetector::HGCEE) layer = HGCEEDetId(rawId()).layer();
	else if (subdetId() == ForwardSubdetector::HGCHEF ||
		 subdetId() == ForwardSubdetector::HGCHEB) {
	  layer = HGCHEDetId(rawId()).layer();
	}
      } 
      setLayer (layer);
    }


  class PileupDataContainer {
  public:
    PileupDataContainer (const vector<double>& etaBinsVec);
    double value (LayerDetId id, double eta) const;
    void add (LayerDetId id, double eta, double value);
    void reset ();
    bool layerAvailable (LayerDetId id) const;
    const vector<double>& getEtaBins () const;
  private:
    vector <double> etaBins;
    map<LayerDetId, vector<double> > container;
    size_t iEta (double eta) const;
  };
  
  PileupDataContainer::PileupDataContainer (const vector<double>& etaBinsVec) 
    : etaBins (etaBinsVec) {}
  
  double PileupDataContainer::value (LayerDetId id, double eta) const {
    auto ptr = container.find (id);
    if (ptr == container.end()) return 0; // nothing available
    return ptr->second [iEta(eta)];
  }
  
  void PileupDataContainer::add (LayerDetId id, double eta, double value) {
    auto ptr = container.find (id);
    if (ptr == container.end()) { // new layer
      container [id] = vector<double> (etaBins.size(), 0.);
    }
    container [id][iEta(eta)] += value;
  }
  
  void PileupDataContainer::reset () {container.clear();}
  
  bool PileupDataContainer::layerAvailable (LayerDetId id) const {
    return container.find (id) != container.end();
  }
  
  const vector<double>& PileupDataContainer::getEtaBins () const {return etaBins;}
  
  size_t PileupDataContainer::iEta (double eta) const {
    auto ptr = lower_bound (etaBins.begin(), etaBins.end(), eta);
    if (ptr == etaBins.end()) return etaBins.size()-1;
    return ptr - etaBins.begin();
  }
}

class PFRecHitPileupSubtractor: public  edm::EDProducer {
public:
  PFRecHitPileupSubtractor (const edm::ParameterSet&);
  virtual ~PFRecHitPileupSubtractor ();
  virtual void produce (edm::Event&, const edm::EventSetup&) override;
  virtual void beginRun (const edm::Run&, const edm::EventSetup&) override;
private:
  edm::InputTag fInputLabel;
  //vector<DetId> fDetectors;
  PileupDataContainer fSquares;
  GenericCaloGeometry fGeometries;
};

PFRecHitPileupSubtractor::PFRecHitPileupSubtractor (const edm::ParameterSet& iConfig)
  : fInputLabel (iConfig.getParameter<edm::InputTag>("input")),
    fSquares (iConfig.getParameter<vector<double> >("etaBins"))
{
  vector<string> detectors (iConfig.getParameter<vector<string> >("detectors"));
  for (auto& detector : detectors) {
    if (detector == "HB") fGeometries.addDetector (DetId (DetId::Hcal, HcalSubdetector::HcalBarrel));
    else if (detector == "HE") fGeometries.addDetector (DetId (DetId::Hcal, HcalSubdetector::HcalEndcap));
    else if (detector == "HF") fGeometries.addDetector (DetId (DetId::Hcal, HcalSubdetector::HcalForward));
    else if (detector == "EB") fGeometries.addDetector (DetId (DetId::Ecal, EcalSubdetector::EcalBarrel));
    else if (detector == "EK") fGeometries.addDetector (DetId (DetId::Ecal, EcalSubdetector::EcalShashlik));
    else if (detector == "HGCEE") fGeometries.addDetector (DetId (DetId::Forward, ForwardSubdetector::HGCEE));
    else if (detector == "HGCHEF") fGeometries.addDetector (DetId (DetId::Forward, ForwardSubdetector::HGCHEF));
    else if (detector == "HGCHEB") fGeometries.addDetector (DetId (DetId::Forward, ForwardSubdetector::HGCHEB));
  }
  produces<RecoPFRecHitRefCandidateCollection>();
}

PFRecHitPileupSubtractor::~PFRecHitPileupSubtractor () {
}

void PFRecHitPileupSubtractor::beginRun (const edm::Run& iRun, const edm::EventSetup& iSetup) {
  fGeometries.reset (iSetup);
  cout << "PFRecHitPileupSubtractor::beginRun detectors: " << fGeometries.allDetectors ().size() << endl;
  for (auto& detector : fGeometries.allDetectors ()) {
    cout << "PFRecHitPileupSubtractor::beginRun detector " << detector<<"->"<<fGeometries.getValidDetIds (detector).size()<<" cells"<< endl;
    for (auto& validDet : fGeometries.getValidDetIds (detector)) {
      //cout << "PFRecHitPileupSubtractor::beginRun 4" << endl;
      fSquares.add (validDet, fGeometries.eta(validDet), fGeometries.etaPhiSquare(validDet));
      if (validDet.det() == DetId::Forward && validDet.subdetId() == ForwardSubdetector::HGCHEF) {
	//cout <<"PFRecHitPileupSubtractor::beginRun-> HGCHEF cell " << HGCHEDetId(validDet.rawId()) << endl;
      }
    }
  }
}

void PFRecHitPileupSubtractor::produce (edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  edm::Handle<RecoPFRecHitRefCandidateCollection> inputHits;
  iEvent.getByLabel(fInputLabel, inputHits);
  //  std::cout << "PFRecHitPileupSubtractor::produce-> input " << inputHits->size() << std::endl;
  std::auto_ptr<RecoPFRecHitRefCandidateCollection> cleanedHits  (new RecoPFRecHitRefCandidateCollection ());
  // 1 loop fill energies
  PileupDataContainer energies (fSquares.getEtaBins ());
  PileupDataContainer energies2 (fSquares.getEtaBins ());
  for (auto& inHit : *inputHits) {
    DetId det (inHit.pfRecHit()->detId());
    LayerDetId ldi (det);
    double eta = fGeometries.eta (det);
    //    double square = fGeometries.etaPhiSquare(inHit.pfRecHit()->detId());
    double p = inHit.p();
    energies.add (ldi, eta, p);
  }
  // 2 loop subtract offsets
  for (auto& inHit : *inputHits) {
    RecoPFRecHitRefCandidate outHit (inHit);
    DetId det (inHit.pfRecHit()->detId());
    LayerDetId ldi (det);
    double eta = fGeometries.eta (det);
    double square = fGeometries.etaPhiSquare(inHit.pfRecHit()->detId());
    double pNew = inHit.p();
    double meanP = square/fSquares.value (ldi, eta)*energies.value (ldi, eta);
    pNew -= meanP;
    if (pNew > 0) {
      double scale = pNew / inHit.p();
      math::PtEtaPhiMLorentzVector newP4 (scale*inHit.pt(),inHit.eta(),inHit.phi(),0);
      outHit.setP4 (newP4);
      cleanedHits->push_back (outHit);
    }
  }
  iEvent.put(cleanedHits);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFRecHitPileupSubtractor);
