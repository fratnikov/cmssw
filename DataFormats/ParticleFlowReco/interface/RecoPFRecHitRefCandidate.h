#ifndef RecoCandidate_RecoPFRecHitRefCandidate_h
#define RecoCandidate_RecoPFRecHitRefCandidate_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

namespace reco {

  class RecoPFRecHitRefCandidate : public LeafCandidate {
  public:
    RecoPFRecHitRefCandidate () 
      : LeafCandidate (),
      etaPhiSquare_ (0)
      {}
    
    RecoPFRecHitRefCandidate (const PFRecHitRef& ref, double mass) 
      : LeafCandidate (0, math::PtEtaPhiMLorentzVector (ref->pt(), ref->eta(), ref->phi(), mass)),
      etaPhiSquare_ (0), ref_ (ref)
    {
      init();
    }
    
    RecoPFRecHitRefCandidate (const PFRecHitRef& ref, Charge q, const Point & vtx = Point( 0, 0, 0 ),
			      int pdgId = 0, int status = 0, bool integerCharge = true ) 
      : LeafCandidate (q, math::PtEtaPhiMLorentzVector (ref->pt(), ref->eta(), ref->phi(), 0), vtx, pdgId, status, integerCharge),
      etaPhiSquare_ (0), ref_ (ref)
      {
	init();
      }
    
    template<typename P4> 
      RecoPFRecHitRefCandidate (const PFRecHitRef& ref, Charge q, const P4 & p4, const Point & vtx = Point( 0, 0, 0 ),
				int pdgId = 0, int status = 0, bool integerCharge = true ) 
      : LeafCandidate (q, p4, vtx, pdgId, status, integerCharge),
      etaPhiSquare_ (0), ref_ (ref)
      {
	init();
      }
    
    virtual ~RecoPFRecHitRefCandidate () {}

    void scaleEnergy (double scale) {
      setP4 (PolarLorentzVector (pt()*scale, eta(), phi(), mass()));
    }
    
    void setEnergy (double e) {
      setP4 (PolarLorentzVector (pt()/energy()*e, eta(), phi(), mass()));
    }
    
    PFRecHitRef const & pfRecHit() const {
      return ref_;
    }   

    double etaPhiSquare () const {return etaPhiSquare_;}

  private:
    double etaPhiSquare_;
    PFRecHitRef ref_;

    void init ();
  };
}

#endif
