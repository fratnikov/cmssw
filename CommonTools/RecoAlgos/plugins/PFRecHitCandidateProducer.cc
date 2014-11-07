
#include "CommonTools/RecoAlgos/src/PFRecHitToRefCandidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef CandidateProducer<
          edm::View<reco::PFRecHit>,
          reco::RecoPFRecHitRefCandidateCollection,
          AnySelector,
          converter::helper::CandConverter<reco::PFRecHit>::type
        > PFRecHitRefCandidateProducer;

DEFINE_FWK_MODULE(PFRecHitRefCandidateProducer);
