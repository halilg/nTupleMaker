import FWCore.ParameterSet.Config as cms

# Default settings. Not meant to be modified
nTupler = cms.EDAnalyzer('nTupleMaker',
   # minTracks = cms.untracked.uint32(0),
   addElec = cms.untracked.bool(False),
   addPhot=cms.untracked.bool(False),
   addMuons =cms.untracked.bool(False),
   addJets=cms.untracked.bool(False),
   addMET=cms.untracked.bool(False),
   addBTag=cms.untracked.bool(False),
   # untracked string fname = "somefile.root"
   labelElec = cms.untracked.string("slimmedElectrons"),
   labelPhot = cms.untracked.string("slimmedPhotons"),
   labelMuons = cms.untracked.string("slimmedMuons"),
   labelJets = cms.untracked.string("slimmedJets"),
   labelMET = cms.untracked.string("slimmedMETs"),
   bTagD = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
   dumpHLT=cms.untracked.vstring("")
)


#demo = cms.EDAnalyzer('nTupleMaker'
#)
