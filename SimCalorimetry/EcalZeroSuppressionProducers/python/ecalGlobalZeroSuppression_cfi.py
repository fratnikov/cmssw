import FWCore.ParameterSet.Config as cms

simEcalGlobalZeroSuppression = cms.EDProducer("EcalZeroSuppressionProducer",
                                              digiProducer = cms.string('simEcalUnsuppressedDigis'),
                                              EBdigiCollection = cms.string(''),
                                              EBZSdigiCollection = cms.string(''),
                                              glbBarrelThreshold = cms.untracked.double(0),
                                              EEdigiCollection = cms.string(''),
                                              EEZSdigiCollection = cms.string(''),
                                              glbEndcapThreshold = cms.untracked.double(0),
                                              EKdigiCollection = cms.string(''),
                                              EKZSdigiCollection = cms.string('ekDigis'),
                                              glbShashlikThreshold = cms.untracked.double(5), #ADC counts of sampleMax > pedestal
)
