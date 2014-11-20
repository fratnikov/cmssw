import FWCore.ParameterSet.Config as cms
pfRecHitPileupSubtractor = cms.EDProducer("PFRecHitPileupSubtractor",
                                          input = cms.InputTag('pfCalibratedRecHitRefsForJets'),
                                          etaBins = cms.vdouble (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                                 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                                 4.2, 4.4, 4.6, 4.8, 5.0),
                                          detectors = cms.vstring ('EB','HB','HF','HGCEE','HGCHEF','HGCHEB')
                                          )
