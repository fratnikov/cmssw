import FWCore.ParameterSet.Config as cms
from RecoLocalCalo.CaloRecCandCreator.caloSuperHitProducer_cfi import *
from RecoLocalCalo.CaloRecCandCreator.caloRecHitCandidates_cff import *
caloSuperHitsHGC = cms.Sequence(caloRecHitsCandidatesHGC +
                                caloSuperHitProducerHGC)

caloSuperHitsShashlik = cms.Sequence(caloRecHitsCandidatesShashlik +
                                caloSuperHitProducerShashlik)
