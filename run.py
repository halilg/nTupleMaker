import FWCore.ParameterSet.Config as cms
import sys

ifname='file:/afs/cern.ch/user/h/halil/hbb/CMSSW_9_4_15/src/UserCode/nTupleMaker/ak8_test.root'
ofname='nTuple_test.root'

if len(sys.argv) > 2: ifname=sys.argv[2]
if len(sys.argv) > 3: ofname=sys.argv[3]

print ifname, '->', ofname

skipEvents=1
maxEvents=1

##########################################################

process = cms.Process("nTupler")

# Default values for module parameters
process.load("UserCode.nTupleMaker.nTupleMaker_cfi")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('nTupleMaker')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( maxEvents ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ifname
    ),
    skipEvents=cms.untracked.uint32(skipEvents)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(ofname)
)

##########################################################
process.nTupler.addMuons=False
process.nTupler.addElec=False
process.nTupler.addPhot=False
process.nTupler.addJets=False
process.nTupler.addMET=False
#process.nTupler.labelMuons = "globalMuons"
process.nTupler.labelFatJets = "ak8CaloJets"

process.p = cms.Path(process.nTupler)
