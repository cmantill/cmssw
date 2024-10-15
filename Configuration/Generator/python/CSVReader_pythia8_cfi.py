import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8CSVReaderGun",

    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(True),

    PGunParameters = cms.PSet(
        Filename = cms.string('GeneratorInterface/Pythia8Interface/test/EtaToGammaAp_ApToEE.csv'),
        ParticleID = cms.vint32(11,13),
        MakeDisplaced = cms.bool(False),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        MinPt = cms.double(25.0),
        MaxPt = cms.double(45.0),
        MinEta = cms.double(-2.4),
        MaxEta = cms.double(2.4),
        MinProdRadius = cms.double(0.0),
        MaxProdRadius = cms.double(10.0),
        NumParticlesPerEvent = cms.int32(5)
        ),

    PythiaParameters = cms.PSet(
        py8ZDecaySettings = cms.vstring(  '23:onMode = off', # turn OFF all Z decays
            '23:onIfAny = 15'  # turn ON Z->tautau
            ),
        parameterSets = cms.vstring(  
            #'py8ZDecaySettings' 
            )
        )
    )
