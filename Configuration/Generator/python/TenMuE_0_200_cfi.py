import FWCore.ParameterSet.Config as cms

# Modified from Configuration/Generator/python/SingleMuPt10_cfi.py
generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        MaxE = cms.double(200.0),
        MinE = cms.double(0.0),
        PartID = cms.vint32(-13,-13,-13,-13,-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359) ## in radians

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('Ten mu e 0 to 200'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)
