import FWCore.ParameterSet.Config as cms

## Default Parameter Sets
#from RecoJets.JetProducers.ak4CaloJets_cfi import ak4CaloJets as _ak4CaloJets

events=10
dpath="/afs/cern.ch/user/h/halil/data/public/Hbb/BulkGravTohhTohbbhbb_width0p10_M-4200_TuneCP2_13TeV-madgraph_pythia8/AODSIM"
rfile="F6AEDBE6-4A53-E911-B821-6C3BE5B59200.root"
ifname="file://"+dpath+"/"+rfile
ofname="ak8_test.root"

#set up a process
processName = "RECO1"
process = cms.Process(processName)

## Calo Towers
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

process.caloTowers = cms.EDProducer("CaloTowerCandidateCreator",
    src = cms.InputTag("towerMaker"),
    e = cms.double(0.0),
    verbose = cms.untracked.int32(0),
    pt = cms.double(0.0),
    minimumE = cms.double(0.0),
    minimumEt = cms.double(0.0),
    et = cms.double(0.0)
)

# jet reconstruction
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#process.load("RecoJets/Configuration/RecoJetsGlobal_cff")

process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("RecoJets.Configuration.CaloTowersRec_cff")
process.load("RecoJets.JetProducers.AnomalousCellParameters_cfi")
process.load("RecoJets.JetProducers.ak4CaloJets_cfi")
process.towerMaker.hbheInput = cms.InputTag("reducedHcalRecHits","hbhereco")
process.towerMaker.hoInput = cms.InputTag("reducedHcalRecHits","horeco")
process.towerMaker.hfInput = cms.InputTag("reducedHcalRecHits","hfreco")
process.towerMaker.ecalInputs = cms.VInputTag(cms.InputTag("reducedEcalRecHitsEE"), cms.InputTag("reducedEcalRecHitsEB"))

process.towerMakerWithHO.hbheInput = cms.InputTag("reducedHcalRecHits","hbhereco")
process.towerMakerWithHO.hoInput = cms.InputTag("reducedHcalRecHits","horeco")
process.towerMakerWithHO.hfInput = cms.InputTag("reducedHcalRecHits","hfreco")
process.towerMakerWithHO.ecalInputs = cms.VInputTag(cms.InputTag("reducedEcalRecHitsEE"), cms.InputTag("reducedEcalRecHitsEB"))

#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

#process.load("Configuration.StandardSequences.Reconstruction_cff")


process.ak8CaloJets=process.ak4CaloJets.clone(rParam = cms.double(0.8))


## Default Sequence
#process.ak8CaloJetMaker = cms.Sequence(
#    process.caloTowersRec*process.caloTowers*process.ak8CaloJets
#    )


process.ak8CaloJetMaker = cms.Sequence(
    process.caloTowersRec*process.caloTowers*process.ak8CaloJets
    )

process.source = cms.Source ("PoolSource", fileNames=cms.untracked.vstring(ifname))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32 (events))

process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string(dpath+ofname),
                outputCommands = cms.untracked.vstring(
                'drop *',
                'keep *_genParticles_*_*',
                'keep *_generator_*_*',
                'keep *_source_*_*',
                'keep *_ak8CaloJets_*_*',
                'keep *_caloTowers_*_*',
                'keep *_towerMaker_*_*',
                'keep *_towerMakerWithHO_*_*',
                'keep *_*_hbhereco_*'
        )
)

# Defines which modules and sequences to run
process.mypath = cms.Path(process.ak8CaloJetMaker)


# A list of analyzers or output modules to be run after all paths have been run.
process.outpath = cms.EndPath(process.out)
