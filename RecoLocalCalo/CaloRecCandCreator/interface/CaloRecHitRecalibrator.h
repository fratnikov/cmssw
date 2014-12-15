#ifndef CaloRecCandCreator_CaloRecHitRecalibrator_h
#define CaloRecCandCreator_CaloRecHitRecalibrator_h

#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/DetId/interface/DetId.h"

class CaloRecHitRecalibrator {
 public:
  virtual double recalibrateEnergy (double energy, double eta, DetId detId) const = 0;
  virtual ~CaloRecHitRecalibrator () {}
  static CaloRecHitRecalibrator* newCaloRecHitRecalibrator (const edm::ParameterSet& pSet);
};

class CaloHGCRecHitRecalibrator : public CaloRecHitRecalibrator {
 public:
  CaloHGCRecHitRecalibrator (const edm::ParameterSet& pSet); 
  virtual ~CaloHGCRecHitRecalibrator () {}
  virtual double recalibrateEnergy (double energy, double eta, DetId detId) const;
 private:
  double mipValueInGeV;
  double coef_a;
  std::vector<double> weights;
  double threshold_MIP;
};

class CaloRecHitScaleRecalibrator : public CaloRecHitRecalibrator {
 public:
  CaloRecHitScaleRecalibrator (const edm::ParameterSet& pSet); 
  virtual ~CaloRecHitScaleRecalibrator () {}
  virtual double recalibrateEnergy (double energy, double eta, DetId detId) const {return energy * scale;}
 private:
  double scale;
};

class CaloRecHitTrivialRecalibrator : public CaloRecHitRecalibrator {
 public:
  CaloRecHitTrivialRecalibrator (const edm::ParameterSet& pSet) {}; 
  virtual ~CaloRecHitTrivialRecalibrator () {}
  virtual double recalibrateEnergy (double energy, double eta, DetId detId) const {return energy;}
};

#endif
