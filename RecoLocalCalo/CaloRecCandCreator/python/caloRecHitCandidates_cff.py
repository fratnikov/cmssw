import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.CaloRecCandCreator.caloRecHitCandidates_cfi import *

caloRecHitsCandidatesHGC = cms.Sequence(caloRecHitCandidatesEB+
                                        caloRecHitCandidatesHBHE+caloRecHitCandidatesHF+
                                        caloRecHitCandidatesHGCEE+caloRecHitCandidatesHGCHEF+caloRecHitCandidatesHGCHEB+
                                        caloRecHitCandidatesHGC
)

caloRecHitsCandidatesShashlik = cms.Sequence(caloRecHitCandidatesEB+
                                        caloRecHitCandidatesHBHE+caloRecHitCandidatesHF+
                                        caloRecHitCandidatesEK+
                                        caloRecHitCandidatesShashlik
)

