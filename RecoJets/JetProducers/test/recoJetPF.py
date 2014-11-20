# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:upgradePLS3 -n 10 --eventcontent FEVTDEBUGHLT,DQM -s RAW2DIGI,L1Reco,RECO,VALIDATION,DQM --datatier GEN-SIM-RECO,DQMIO --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco --magField 38T_PostLS1 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('USER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_6_2_0_SLHC20/RelValSinglePiE50HCAL/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV5-v1/00000/0A797444-AB5F-E411-98ED-002618943885.root',
'/store/relval/CMSSW_6_2_0_SLHC20/RelValSinglePiE50HCAL/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV5-v1/00000/BE5E324C-AB5F-E411-8740-0025905964B4.root')
#'/store/relval/CMSSW_6_2_0_SLHC20/RelValSingleGammaPt35Extended/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV5-v1/00000/DEBEC36D-A05F-E411-BA3E-002618943865.root',
#'/store/relval/CMSSW_6_2_0_SLHC20/RelValSingleGammaPt35Extended/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV5-v1/00000/EC93226D-A05F-E411-A5C9-002618943940.root')
#    fileNames = cms.untracked.vstring('file:/uscms/home/ratnikov/cms/SLHC/jets/CMSSW_6_2_0_SLHC19/src/step3.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    ##outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    outputCommands = cms.untracked.vstring("keep *_*_*_*"),                                         
    fileName = cms.untracked.string('file:data/RelValPiE50.root'),
#    fileName = cms.untracked.string('file:data/RelValSingleGammaPt35.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

#process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.DQMEventContent.outputCommands,
#    fileName = cms.untracked.string('file:step3_inDQM.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('DQMIO')
#    )
#)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.prevalidation_step = cms.Path(process.prevalidation)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.validation_step = cms.EndPath(process.validation)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)

#------- new stuff -------
process.pfRecHitPileupSubtractor = cms.EDProducer("PFRecHitPileupSubtractor",
                                          input = cms.InputTag('pfCalibratedRecHitRefsForJets'),
                                          etaBins = cms.vdouble (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                                 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                                 4.2, 4.4, 4.6, 4.8, 5.0, 999.),
#                                          detectors = cms.vstring ('HE')
                                          detectors = cms.vstring ('HGCEE', 'HGCHEF', 'HGCHEB')
#                                          detectors = cms.vstring ('HGCHEF')
)
process.pileupSubtractor = cms.Sequence(process.pfRecHitPileupSubtractor)
process.pileupSubtractor_step = cms.Path (process.pileupSubtractor)

process.load ('RecoJets.Configuration.RecoPFRecHitJets_cff')
process.jetreco_step = cms.Path (process.particleFlowCluster*process.recoPFRecHitJets)
#process.schedule = cms.Schedule(process.pileupSubtractor_step, process.jetreco_step, process.FEVTDEBUGHLToutput_step)
process.schedule = cms.Schedule(process.jetreco_step, process.FEVTDEBUGHLToutput_step)


#------- new stuff end -------

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.prevalidation_step,process.validation_step,process.dqmoffline_step,process.FEVTDEBUGHLToutput_step,process.DQMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon 

#call to customisation function cust_2023HGCalMuon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCalMuon(process)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
