#ifndef RecoAlgos_PFClusterToRefCandidate_h
#define RecoAlgos_PFClusterToRefCandidate_h
#include "CommonTools/RecoAlgos/src/MassiveCandidateConverter.h"
#include "CommonTools/RecoAlgos/src/CandidateProducer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

namespace converter {

  struct PFRecHitToRefCandidate : public MassiveCandidateConverter {
    typedef reco::PFRecHit value_type;
    typedef reco::PFRecHitCollection Components;
    typedef reco::RecoPFRecHitRefCandidate Candidate;
    PFRecHitToRefCandidate(const edm::ParameterSet & cfg) : 
      MassiveCandidateConverter(cfg) {
    }
    void convert(reco::PFRecHitRef pfclusterRef, reco::RecoPFRecHitRefCandidate & c) const {
      c = reco::RecoPFRecHitRefCandidate( pfclusterRef, sqrt(massSqr_) );
    }  
  };

  namespace helper {
    template<>
    struct CandConverter<reco::PFRecHit> { 
      typedef PFRecHitToRefCandidate type;
    };
  }

}

#endif
