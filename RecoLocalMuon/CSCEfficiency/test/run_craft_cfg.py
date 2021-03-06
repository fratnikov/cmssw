import FWCore.ParameterSet.Config as cms


process = cms.Process("PROC")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'CRAFT_30X::All'
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
# for Beam
#process.load("Configuration.StandardSequences.Reconstruction_cff")
# for Cosmics
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")


from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.source = cms.Source ("PoolSource",
fileNames = cms.untracked.vstring (
       '/store/data/Commissioning08/Cosmics/RECO/CRAFT_ALL_V4_ReReco-v1/0026/321C18E5-86C6-DD11-82FE-0019B9E4905D.root'
),
#skipEvents = cms.untracked.uint32(808)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.ana = cms.EDFilter("CSCEfficiency",
    MuonServiceProxy,
#    printout_NEvents = cms.untracked.uint32(1000),
#    rootFileName = cms.untracked.string('cscHists.root'),
#    getAbsoluteEfficiency = cms.untracked.bool(True),    
    useDigis = cms.untracked.bool(True),
    runOnData = cms.untracked.bool(True),
#    IPdata = cms.untracked.bool(False),
#    Beamdata = cms.untracked.bool(False),
#    distanceFromDeadZone = cms.untracked.double(10.0),
    applyIPangleCuts = cms.untracked.bool(True),
    alctDigiTag = cms.InputTag("muonCSCDigis","MuonCSCALCTDigi"),
    clctDigiTag = cms.InputTag("muonCSCDigis","MuonCSCCLCTDigi"),
    corrlctDigiTag = cms.InputTag("muonCSCDigis","MuonCSCCorrelatedLCTDigi"),
    stripDigiTag = cms.InputTag("muonCSCDigis","MuonCSCStripDigi"),
    wireDigiTag = cms.InputTag("muonCSCDigis","MuonCSCWireDigi"),
    rechitDigiTag = cms.InputTag("csc2DRecHits"),
    segmentDigiTag = cms.InputTag("cscSegments"),
    simHitTag = cms.InputTag("g4SimHits", "MuonCSCHits"),
#    tracksTag = cms.InputTag("globalMuons"),
#("generalTracks"), ("standAloneMuons"), ("standAloneMuons","UpdatedAtVtx"), ("globalMuons")
#    tracksTag = cms.InputTag("generalTracks"),
#cosmics: cosmicMuonsEndCapsOnly, cosmicMuons, globalCosmicMuons, ctfWithMaterialTracksP5
    tracksTag = cms.InputTag("cosmicMuons"),
#    maxNormChi2 = cms.untracked.double(3.0),
#    minTrackHits = cms.untracked.uint32(10),
# if no magnetic filed - P values don't matter
#    minP = cms.untracked.double(10.0),
#    maxP = cms.untracked.double(100.0)
# trigger
    useTrigger = cms.untracked.bool(True),
    HLTriggerResults = cms.InputTag('TriggerResults','','HLT' ),
    myTriggers = cms.vstring("HLT_L1MuOpen"),#, "HLT_L1_CSCMuonHalo"),
#    myTriggers = cms.vstring("HLT_L1_CSCMuonHalo"),
    andOr = cms.untracked.bool(False), # "true" means OR 

)

# memory issues when too many runs - a cure: 
#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*","drop *_lumiProducer_*_REPACKER")

# for running Iguana
process.VisConfigurationService = cms.Service("VisConfigurationService",
    EnabledTwigs = cms.untracked.vstring('/Objects/CMS Event and Detector/Muon/Endcap/CSCs','*cosmicMuons*'),
    VisAutoStart = cms.untracked.bool(False),
    ContentProxies = cms.untracked.vstring('Reco/Calorimetry',
        'Reco/Candidate',
        'Reco/Detector',
#        'Reco/Tracker', 
        'Reco/Muon',
        'Reco/MuonCSC',
        'Reco/MuonCSC/Strip digis',
        'Reco/MuonCSC/Wire digis',
        'Reco/MuonCSC/Rec Hit 2D',
        'Reco/MuonCSC/Segments',
        'Reco/MuonDT',
        'Reco/MuonRPC',
        'Reco/Tools',
#        'Reco/Trigger',
        'Simulation/Hits',
                                   )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/uscms_data/d1/stoyan/data/tmp/test_sample.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('analyze')
    )
)
# Uncomment for output skim (and change fileName path above)
#process.outpath = cms.EndPath(
#  process.out 
#)

process.analyze = cms.Path(process.ana)
#process.p = cms.Path(process.muonCSCDigis*process.csc2DRecHits*process.cscSegments*process.ana)

    
