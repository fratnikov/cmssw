#include "DataFormats/RecoCandidate/interface/CaloRecHitCandidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/Merger.h"

typedef Merger<std::vector<reco::CaloRecHitCandidate> > CaloRecHitCandidateMerger;

DEFINE_FWK_MODULE(CaloRecHitCandidateMerger);
