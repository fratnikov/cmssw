#include "DataFormats/ParticleFlowReco/interface/RecoPFRecHitRefCandidateFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/Merger.h"

typedef Merger<reco::RecoPFRecHitRefCandidateCollection>
        PFRecHitRefCandidateMerger;

DEFINE_FWK_MODULE(PFRecHitRefCandidateMerger);
