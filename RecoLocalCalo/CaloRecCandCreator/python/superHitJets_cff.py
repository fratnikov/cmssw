import FWCore.ParameterSet.Config as cms

#from RecoJets.JetProducers.ak5PFRecHitJets_cfi import ak5PFRecHitJets
##from RecoJets.JetProducers.PFClusterJetParameters_cfi import *

SuperHitJetParameters = cms.PSet( 
    src            = cms.InputTag('caloSuperHitProducerHGC'),
    srcPVs         = cms.InputTag('offlinePrimaryVertices'),
    jetType        = cms.string('BasicJet'),
    doOutputJets   = cms.bool(True),
    # minimum jet pt
    jetPtMin       = cms.double(3.0),
    # minimum calo tower input et
#    inputEtMin     = cms.double(0.3),
    inputEtMin     = cms.double(0.00001),
    # minimum calo tower input energy
    inputEMin      = cms.double(0.0),
    # primary vertex correction
    doPVCorrection = cms.bool(True),
    # pileup with offset correction
    doPUOffsetCorr = cms.bool(False),
    # if pileup is false, these are not read:
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),  
    # fastjet-style pileup 
    doAreaFastjet       = cms.bool( False),
    doRhoFastjet        = cms.bool( False),
    doAreaDiskApprox    = cms.bool( False),
    Active_Area_Repeats = cms.int32(    1),
    GhostArea           = cms.double(0.01),
    Ghost_EtaMax        = cms.double( 5.0),
    Rho_EtaMax          = cms.double( 4.4),
    voronoiRfact        = cms.double(-0.9),
    useDeterministicSeed= cms.bool( True ),
    minSeed             = cms.uint32( 14327 )
)


from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

ak5SuperHitJets = cms.EDProducer(
    "FastjetJetProducer",
    SuperHitJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.5)
    )

ak5RecHitJets = cms.EDProducer(
    "FastjetJetProducer",
    SuperHitJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.5)
    )
ak5RecHitJets.src = cms.InputTag('caloRecHitCandidatesHGC')


from RecoLocalCalo.CaloRecCandCreator.caloSuperHitProducer_cff import *

recoCaloJets = cms.Sequence(caloSuperHitsHGC+ak5SuperHitJets+ak5RecHitJets)

