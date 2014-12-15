#include "RecoLocalCalo/CaloRecCandCreator/interface/CaloRecHitRecalibrator.h"

#include <iostream>

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

using namespace std;

CaloRecHitRecalibrator* CaloRecHitRecalibrator::newCaloRecHitRecalibrator (const edm::ParameterSet& pSet) {
  string label (pSet.getParameter<string>("label"));
  if (label == "HGCalRecalibrator") return new CaloHGCRecHitRecalibrator (pSet);
  else if (label == "ScaleRecalibrator") return new CaloRecHitScaleRecalibrator (pSet);
  else return new CaloRecHitTrivialRecalibrator (pSet);
}

CaloHGCRecHitRecalibrator:: CaloHGCRecHitRecalibrator (const edm::ParameterSet& pSet) 
  : mipValueInGeV (pSet.getParameter<double>("MipValueInGeV")),
    coef_a (pSet.getParameter<double>("effMip_to_InverseGeV_a")),	
    weights (pSet.getParameter<vector<double> >("weights")),
    threshold_MIP (pSet.getParameter<double>("thresholdInMIPs"))
{ 
}

double  CaloHGCRecHitRecalibrator::recalibrateEnergy (double energy, double eta, DetId detId) const {
  if (DetId (detId).det() == DetId::Forward) {
    int subdetId = detId.subdetId();
    int layer = subdetId == ForwardSubdetector::HGCEE ? HGCEEDetId(detId).layer() :
      (subdetId == ForwardSubdetector::HGCHEF || ForwardSubdetector::HGCHEB) ? HGCHEDetId(detId).layer() :
      -1;
    if (layer >= 1 || layer <= int(weights.size())) {
      double energy_MIP = energy/mipValueInGeV;
      double etaCorrection = fabs(tanh(eta)); // |pz/p|
      if (energy_MIP < threshold_MIP / etaCorrection) return 0; // under the threshold
      double eCorr = weights[layer]*energy_MIP;
      //double effMIP_to_InvGeV = coef_a/(1.0 + std::exp(-coef_c - coef_b*std::cosh(eta)));
      double effMIP_to_InvGeV = coef_a;
      return eCorr/effMIP_to_InvGeV;
    }
    cout << "HGCalRecalibrator::recalibrateEnergy-> invalig layer " << layer << " expect 1.." << weights.size() << endl;
  }
  else {
    cout << "HGCalRecalibrator::recalibrateEnergy-> invalig detector:  Forward is expected" << endl;
  }
  return 0;
}

CaloRecHitScaleRecalibrator::CaloRecHitScaleRecalibrator (const edm::ParameterSet& pSet) :
  scale (pSet.getParameter<double>("scale"))
{} 





