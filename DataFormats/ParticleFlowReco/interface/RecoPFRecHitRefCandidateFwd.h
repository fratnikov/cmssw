#ifndef RecoCandidate_RecoPFRecHitRefCandidateFwd_h
#define RecoCandidate_RecoPFRecHitRefCandidateFwd_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidate.h"

namespace reco {

  /// collectin of LeafRefCandidateT<reco::TrackRef>  objects
  typedef std::vector<RecoPFRecHitRefCandidate > RecoPFRecHitRefCandidateCollection;

  /// reference to an object in a collection of RecoPFRecHitRefCandidate objects
  typedef edm::Ref<RecoPFRecHitRefCandidateCollection> RecoPFRecHitRefCandidateRef;

  /// reference to a collection of RecoPFRecHitRefCandidate objects
  typedef edm::RefProd<RecoPFRecHitRefCandidateCollection> RecoPFRecHitRefCandidateRefProd;

  /// vector of objects in the same collection of RecoPFRecHitRefCandidate objects
  typedef edm::RefVector<RecoPFRecHitRefCandidateCollection> RecoPFRecHitRefCandidateRefVector;

  /// iterator over a vector of reference to RecoPFRecHitRefCandidate objects
  typedef RecoPFRecHitRefCandidateRefVector::iterator recoPFRecHitRefCandidate_iterator;  

  typedef edm::RefToBase<reco::Candidate> RecoPFRecHitRefCandidateRefToBase;

/*   /// this needs to go here, it's a class template in the DF/Candidate package */
/*   /// that requires the knowledge of the DF/TrackReco dictionaries */
/*   typedef std::vector<RecoPFRecHitRefCandidateBase> RecoPFRecHitRefCandidateBaseCollection; */
/*   typedef edm::Ref<RecoPFRecHitRefCandidateBaseCollection> RecoPFRecHitRefCandidateBaseRef; */
/*   typedef edm::RefVector<RecoPFRecHitRefCandidateBaseCollection> RecoPFRecHitRefCandidateBaseRefVector; */
/*   typedef edm::RefProd<RecoPFRecHitRefCandidateBaseCollection> RecoPFRecHitRefCandidateBaseRefProd; */
/*   typedef edm::RefToBase<reco::Candidate> RecoPFRecHitRefCandidateBaseRefToBase; */

}

#endif
