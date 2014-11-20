import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

pfRecHitRefsForJetsHGEE = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitHGCEE'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJetsHGHEF = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitHGCHEF'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJetsHGHEB = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitHGCHEB'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJetsHB = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitHBHE'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJetsHF = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitHF'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJetsEB = cms.EDProducer("PFRecHitRefCandidateProducer",
    src          = cms.InputTag('particleFlowRecHitECAL'),
    particleType = cms.string('gamma')
)

pfRecHitRefsForJets = cms.EDProducer("PFRecHitRefCandidateMerger",
    src = cms.VInputTag("pfRecHitRefsForJetsHGEE", "pfRecHitRefsForJetsHGHEF", "pfRecHitRefsForJetsHGHEB",
                        "pfRecHitRefsForJetsHB", "pfRecHitRefsForJetsHF", "pfRecHitRefsForJetsEB")
)

pfCalibratedRecHitRefsForJets = cms.EDProducer("PFRecHitCalibrator",
                                               input = cms.InputTag('pfRecHitRefsForJets'),
                                               HGCEEElectronEnergyCalibrator = cms.PSet (
        #weights for layers from P.Silva (24 October 2014)
        ## this is for V5!
        # MIP effective to 1.0/GeV (from fit to data of P. Silva)
        #f(x) = a/(1-exp(-bx - c))
        # x = cosh(eta)
        # a = 82.8
        # b = 1e6
        # c = 1e6
        weights = cms.vdouble([0]+[0.080]+[0.62]*9+[0.81]*9+[1.19]*8),
        effMip_to_InverseGeV_a = cms.double(82.8),
        effMip_to_InverseGeV_b = cms.double(1e6),
        effMip_to_InverseGeV_c = cms.double(1e6),
        MipValueInGeV = cms.double(55.1*1e-6)
        ),
                                               HGCHEFHadronicEnergyCalibrator = cms.PSet (
        weights = cms.vdouble([0]+[0.0464]+[0.0474]*11),
        effMip_to_InverseGeV_a = cms.double(1.0),
        effMip_to_InverseGeV_b = cms.double(1e6),
        effMip_to_InverseGeV_c = cms.double(1e6),
        MipValueInGeV = cms.double(85.0*1e-6)
        ),
                                               HGCHEBHadronicEnergyCalibrator = cms.PSet (
        weights = cms.vdouble([0]+[0.1215]*12),
        effMip_to_InverseGeV_a = cms.double(1.0),
        effMip_to_InverseGeV_b = cms.double(1e6),
        effMip_to_InverseGeV_c = cms.double(1e6),
        MipValueInGeV = cms.double(1498.4*1e-6)
        )
                                               
                                               )
