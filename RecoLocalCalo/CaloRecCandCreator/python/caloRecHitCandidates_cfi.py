import FWCore.ParameterSet.Config as cms

caloRecHitCandidatesEB = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('None')
           )
)

caloRecHitCandidatesHBHE = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag('hbheUpgradeReco'),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('None')
           )
)

caloRecHitCandidatesHF = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag('hfUpgradeReco'),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('None')
           )
)

caloRecHitCandidatesEK = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag("ecalRecHit","EcalRecHitsEK"),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('None')
           )
)

caloRecHitCandidatesHGCEE = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag("HGCalRecHit","HGCEERecHits"),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('HGCalRecalibrator'),
              weights = cms.vdouble([0]+[0.080]+[0.62]*10+[0.809]*10+[1.239]*9),
              effMip_to_InverseGeV_a = cms.double(82.8),
              MipValueInGeV = cms.double(55.1*1e-6)
           )
)

caloRecHitCandidatesHGCHEF = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('HGCalRecalibrator'),
              weights = cms.vdouble([0]+[0.0464]+[0.0474]*11),
              effMip_to_InverseGeV_a = cms.double(1.0),
              MipValueInGeV = cms.double(85.0*1e-6)
           )
)

caloRecHitCandidatesHGCHEB = cms.EDProducer("CaloRecHitCandidateProducer",
           input = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
           energyThreshold = cms.double(0),
           recalibrator = cms.PSet (
              label = cms.string('HGCalRecalibrator'),
              weights = cms.vdouble([0]+[0.1215]*12),
              effMip_to_InverseGeV_a = cms.double(1.0),
              MipValueInGeV = cms.double(1498.4*1e-6)
           )
)

caloRecHitCandidatesHGC = cms.EDProducer("CaloRecHitCandidateMerger",
              src = cms.VInputTag("caloRecHitCandidatesEB",
                                  "caloRecHitCandidatesHBHE",
                                  "caloRecHitCandidatesHF",
                                  "caloRecHitCandidatesHGCEE",
                                  "caloRecHitCandidatesHGCHEF",
                                  "caloRecHitCandidatesHGCHEB"
                                  )
)
                                  
caloRecHitCandidatesShashlik = cms.EDProducer("CaloRecHitCandidateMerger",
              src = cms.VInputTag("caloRecHitCandidatesEB",
                                  "caloRecHitCandidatesHBHE",
                                  "caloRecHitCandidatesHF",
                                  "caloRecHitCandidatesEK",
                                  )
)
                                  

