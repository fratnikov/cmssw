#ifndef RecoCandidate_CaloSuperHitFwd_h
#define RecoCandidate_CaloSuperHitFwd_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/ParticleFlowReco/interface/CaloSuperHit.h"

namespace reco {

  /// collectin of LeafRefCandidateT<reco::TrackRef>  objects
  typedef std::vector<CaloSuperHit > CaloSuperHitCollection;

  /// reference to an object in a collection of CaloSuperHit objects
  typedef edm::Ref<CaloSuperHitCollection> CaloSuperHitRef;

  /// reference to a collection of CaloSuperHit objects
  typedef edm::RefProd<CaloSuperHitCollection> CaloSuperHitRefProd;

  /// vector of objects in the same collection of CaloSuperHit objects
  typedef edm::RefVector<CaloSuperHitCollection> CaloSuperHitRefVector;

  /// iterator over a vector of reference to CaloSuperHit objects
  typedef CaloSuperHitRefVector::iterator CaloSuperHit_iterator;  

  typedef edm::RefToBase<reco::Candidate> CaloSuperHitRefToBase;
}

#endif
