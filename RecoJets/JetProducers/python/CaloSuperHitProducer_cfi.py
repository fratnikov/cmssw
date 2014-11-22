import FWCore.ParameterSet.Config as cms
caloSuperHitProducer = cms.EDProducer("CaloSuperHitProducer",
                                          input = cms.InputTag('pfCalibratedRecHitRefsForJets'),
                                          etaBins = cms.vdouble ([(i+1)*0.087 for i in range(int(1.7/0.087+0.5))]+
                                                                 [(i+1)*0.087 for i in range(int(1.7/0.087+0.5), int(3.0/0.087+0.5))]+
                                                                 [(i+1)*0.087 for i in range(int(3.0/0.087+0.5), int(5.0/0.087+0.5))]
                                                                 ),
                                          phiBins = cms.vuint32 ([72 for i in range(int(1.7/0.087+0.5))]+
                                                                 [144 for i in range(int(1.7/0.087+0.5), int(3.0/0.087+0.5))]+
                                                                 [36 for i in range(int(3.0/0.087+0.5), int(5.0/0.087+0.5))]
                                                                 ),
                                          detectors = cms.vstring ('EB','HB','HF','HGCEE','HGCHEF','HGCHEB')
                                          )
