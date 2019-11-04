import os
import FWCore.ParameterSet.Config as cms

####### Configuration #######
processName = "RECO1"
events=10
ifpath="/afs/cern.ch/user/h/halil/data/public/Hbb/BulkGravTohhTohbbhbb_width0p10_M-4200_TuneCP2_13TeV-madgraph_pythia8/AODSIM"
ifile="F6AEDBE6-4A53-E911-B821-6C3BE5B59200.root"
ifname="file://"+ os.path.join(ifpath, ifile)
ofpath=""
ofname="ak8_test.root"
outputcmd = (   'drop *',
                'keep *_genParticles_*_*',
                'keep *_generator_*_*',
                'keep *_source_*_*',
                'keep *_ak8CaloJets_*_*',
                'keep *_caloTowers_*_*',
                'keep *_towerMaker_*_*',
                'keep *_towerMakerWithHO_*_*',
                'keep *_reducedHcalRecHits_*_*', # temporary
                'keep *_*_hbhereco_*'
            )
####### /Configuration #######

#set up the process
process = cms.Process(processName)
process.source = cms.Source ("PoolSource", fileNames=cms.untracked.vstring(ifname))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32 (events))
process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string(os.path.join(ofpath,ofname)),
                outputCommands = cms.untracked.vstring(*outputcmd)
              )
# commenting out because don't know if actually needed. jet reconstruction seems fine without it
# Calo Tower candidate producer
# process.caloTowers = cms.EDProducer("CaloTowerCandidateCreator",
#     src = cms.InputTag("towerMaker"),
#     e = cms.double(0.0),
#     verbose = cms.untracked.int32(0),
#     pt = cms.double(0.0),
#     minimumE = cms.double(0.0),
#     minimumEt = cms.double(0.0),
#     et = cms.double(0.0)
# )

# General reconstruction
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# CaloTower reconstruction
process.load("RecoJets.Configuration.CaloTowersRec_cff")          # 

# Without these, towermaker cant find the products in the input file. Don't know why.
# This is the error message:
# -----------------------
# Begin processing the 1st record. Run 1, Event 14005, LumiSection 15 at 01-Nov-2019 09:14:54.217 CET
# ----- Begin Fatal Exception 01-Nov-2019 09:14:56 CET-----------------------
# An exception of category 'ProductNotFound' occurred while
#    [0] Processing  Event run: 1 lumi: 15 event: 14005 stream: 0
#    [1] Running path 'mypath'
#    [2] Calling method for module CaloTowersCreator/'towerMaker'
# Exception Message:
# Principal::getByToken: Found zero products matching all criteria
# Looking for type: edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >
# Looking for module label: hbhereco
# Looking for productInstanceName: 
# -----------------------

process.towerMaker.hbheInput = cms.InputTag("reducedHcalRecHits","hbhereco")
process.towerMaker.hoInput = cms.InputTag("reducedHcalRecHits","horeco")
process.towerMaker.hfInput = cms.InputTag("reducedHcalRecHits","hfreco")
process.towerMaker.ecalInputs = cms.VInputTag(cms.InputTag("reducedEcalRecHitsEE"), cms.InputTag("reducedEcalRecHitsEB"))

process.towerMakerWithHO.hbheInput = cms.InputTag("reducedHcalRecHits","hbhereco")
process.towerMakerWithHO.hoInput = cms.InputTag("reducedHcalRecHits","horeco")
process.towerMakerWithHO.hfInput = cms.InputTag("reducedHcalRecHits","hfreco")
process.towerMakerWithHO.ecalInputs = cms.VInputTag(cms.InputTag("reducedEcalRecHitsEE"), cms.InputTag("reducedEcalRecHitsEB"))

# Jet reconstruction
process.load("RecoJets.JetProducers.ak4CaloJets_cfi")
process.ak8CaloJets=process.ak4CaloJets.clone(rParam = cms.double(0.8))

# Put caloTower reco and jet reco into a single sequence
process.mySeq = cms.Sequence(
    process.caloTowersRec * 
    #process.caloTowers * # this makes the product: edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >    "caloTowers"
    process.ak8CaloJets
    )

# Final steps
process.mypath = cms.Path(process.mySeq)
process.outpath = cms.EndPath(process.out)
