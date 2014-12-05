#ifndef RecoCandidate_CaloSuperHit_h
#define RecoCandidate_CaloSuperHit_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace reco {
  
  class CaloSuperHit : public LeafCandidate {
  public:
    CaloSuperHit () 
      : LeafCandidate (),
      etaCenter_ (0),
      phiCenter_ (0),
      dEta_(0),
      dPhi_(0)
	{}
    
    CaloSuperHit (double etaCenter, double phiCenter, double dEta, double dPhi)
      : LeafCandidate (),
      etaCenter_ (etaCenter),
      phiCenter_ (phiCenter),
      dEta_(dEta),
      dPhi_(dPhi)
      {}
    
    virtual ~CaloSuperHit () {}
      
    template<typename P4> 
      void addHit (DetId idAdd, const P4& p4Add) {
      constituents_.push_back (idAdd);
      LeafCandidate::LorentzVector pSum = p4Add;
      pSum += p4();
      setP4 (pSum);
      setMass (0);
    }
    
    bool belongs (double eta, double phi) const {
      if (deltaPhi (phiCenter (), phi) > dPhi()) return false;
      if (fabs(etaCenter () - eta) > dEta()) return false;
      return true;
    }

    double etaPhiSquare () const {
      return 4*dEta ()*dPhi ();
    }

    const std::vector<DetId>& constituents () const {return constituents_;}
    double etaCenter () const {return etaCenter_;}
    double phiCenter () const {return phiCenter_;}
    double dEta () const {return dEta_;}
    double dPhi () const {return dPhi_;}

  private:
    double etaCenter_;
    double phiCenter_;
    double dEta_;
    double dPhi_;
    std::vector<DetId> constituents_;
  };
}

#endif
