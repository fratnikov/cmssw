import FWCore.ParameterSet.Config as cms

shashlikAnalyzer = cms.EDAnalyzer('ShashlikAnalyzer',
#                                  EKdigiCollection = cms.InputTag("simEcalUnsuppressedDigis"),
                                  #EKdigiCollection = cms.InputTag("simEcalDigis","ekDigis"),
                                  EKdigiCollection = cms.InputTag("simEcalGlobalZeroSuppression","ekDigis"),
                                  EKrecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEK"),
)
